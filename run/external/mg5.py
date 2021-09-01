import subprocess
import json
import re
import shutil
import io
import os

import numpy as np
import pandas as pd
import pygit2
import lz4.frame

from .. import paths, tools
from . import interface


class Mg5(interface.External):
    pass


def run_mc(name):
    source = paths.runcards / name
    dest = tools.create_folder(name)
    mg5_dir = dest / name

    # copy the output file to the directory and replace the variables
    output = (source / "output.txt").read_text().replace("@OUTPUT@", name)
    output_file = dest / "output.txt"
    output_file.write_text(output)

    # create output folder
    output_log = tools.run_subprocess([str(paths.mg5_exe), str(output_file)], dest=dest)
    (dest / "output.log").write_text(output_log)

    # copy patches if there are any; use xargs to properly signal failures
    for p in source.iterdir():
        if p.suffix == ".patch":
            subprocess.run(
                "patch -p1".split(), input=p.read_text(), text=True, cwd=mg5_dir
            )

    # enforce proper analysis
    # - copy analysis.f
    analysis = (paths.runcards / name / "analysis.f").read_text()
    (mg5_dir / "FixedOrderAnalysis" / f"{name}.f").write_text(analysis)
    # - update analysis card
    analysis_card = mg5_dir / "Cards" / "FO_analyse_card.dat"
    analysis_card.write_text(
        analysis_card.read_text().replace("analysis_HwU_template", name)
    )

    # copy the launch file to the directory and replace the variables
    launch = (source / "launch.txt").read_text().replace("@OUTPUT@", name)

    # TODO: write a list with variables that should be replaced in the launch file; for the time
    # being we create the file here, but in the future it should be read from the theory database
    variables = json.loads((paths.pkg / "variables.json").read_text())

    # replace the variables with their values
    for name, value in variables.items():
        launch = launch.replace(f"@{name}@", value)

    # perform simple arithmetic on lines containing 'set' and '=' and arithmetic operators
    #  interpreter = asteval.Interpreter()  # use asteval for safety
    #  lines = []
    #  pattern = re.compile(r"(set [\w_]* = )(.*)")
    #  for line in launch.splitlines():
    #  m = re.fullmatch(pattern, line)
    #  if m is not None:
    #  line = m[1] + str(interpreter.eval(m[2]))
    #  lines.append(line)
    #  launch = "\n".join(lines)

    # finally write launch
    launch_file = dest / "launch.txt"
    launch_file.write_text(launch)

    # parse launch file for user-defined cuts
    user_cuts_pattern = re.compile(
        r"^#user_defined_cut set (\w+)\s+=\s+([+-]?\d+(?:\.\d+)?|True|False)$"
    )
    user_cuts = []
    for line in launch.splitlines():
        m = re.fullmatch(user_cuts_pattern, line)
        if m is not None:
            user_cuts.append((m[1], m[2]))

    # if there are user-defined cuts, implement them
    apply_user_cuts(mg5_dir / "SubProcesses" / "cuts.f", user_cuts)

    # parse launch file for user-defined minimum tau
    user_taumin_pattern = re.compile(r"^#user_defined_tau_min (.*)")
    user_taumin = None
    for line in launch.splitlines():
        m = re.fullmatch(user_taumin_pattern, line)
        if m is not None:
            try:
                user_taumin = float(m[1])
            except ValueError:
                raise ValueError("User defined tau_min is expected to be a number")

    if user_taumin is not None:
        set_tau_min_patch = (
            (paths.patches / "set_tau_min.patch")
            .read_text()
            .replace("@TAU_MIN@", f"{user_taumin}d0")
        )
        (dest / "set_tau_min.patch").write_text(set_tau_min_patch)
        subprocess.run(
            f"patch -p1 -d '{mg5_dir}'".split(), input=set_tau_min_patch, text=True
        )

    # parse launch file for other patches
    enable_patches_pattern = re.compile(r"^#enable_patch (.*)")
    enable_patches_list = []
    for line in launch.splitlines():
        m = re.fullmatch(user_taumin_pattern, line)
        if m is not None:
            enable_patches_list.append(m[1])

    if len(enable_patches_list) != 0:
        for patch in enable_patches_list:
            patch_file = paths.patches / patch
            patch_file = patch_file.with_suffix(patch_file.suffix + ".patch")
            if not patch_file.exists():
                raise ValueError(
                    f"Patch '{patch}' requested, but does not exist in patches folder"
                )
            subprocess.run(
                f"patch -p1 -d '{mg5_dir}'".split(),
                input=patch_file.read_text(),
                text=True,
            )

    # launch run
    launch_log = tools.run_subprocess([str(paths.mg5_exe), str(launch_file)], dest=dest)
    (dest / "launch.log").write_text(launch_log)

    return dest


def find_marker_position(insertion_marker, contents):
    marker_pos = -1

    for lineno, value in enumerate(contents):
        if insertion_marker in value:
            marker_pos = lineno
            break

    if marker_pos == -1:
        raise ValueError(
            "Error: could not find insertion marker `{insertion_marker}` in cut file `{file_path}`"
        )

    return marker_pos


def apply_user_cuts(file_path, user_cuts):
    with open(file_path, "r") as fd:
        contents = fd.readlines()

    # insert variable declaration
    marker_pos = find_marker_position("logical function passcuts_user", contents)
    marker_pos = marker_pos + 8

    for fname in paths.cuts_variables.iterdir():
        name = fname.stem
        if any(i[0].startswith(name) for i in user_cuts):
            contents.insert(marker_pos, fname.read_text())

    marker_pos = find_marker_position("USER-DEFINED CUTS", contents)
    # skip some lines with comments
    marker_pos = marker_pos + 4
    # insert and empty line
    contents.insert(marker_pos - 1, "\n")

    for name, value in reversed(user_cuts):
        # map to fortran syntax
        if value == "True":
            value = ".true."
        elif value == "False":
            value = ".false."
        else:
            try:
                float(value)
            except ValueError:
                raise ValueError(f"Error: format of value `{value}` not understood")

            value = value + "d0"

        code = (paths.cuts_code / f"{name}.f").read_text().format(value)
        contents.insert(marker_pos, code)

    with open(file_path, "w") as fd:
        fd.writelines(contents)


def merge(name, dest):
    source = paths.runcards / name
    mg5_dir = dest / name
    grid = dest / f"{name}.pineappl"
    gridtmp = dest / f"{name}.pineappl.tmp"
    pineappl = paths.pineappl_exe()

    # merge the final bins
    mg5_grids = " ".join(
        sorted(str(p) for p in mg5_dir.glob("Events/run_01*/amcblast_obs_*.pineappl"))
    )
    subprocess.run(f"{pineappl} merge {grid} {mg5_grids}".split())

    # optimize the grids
    subprocess.run(f"{pineappl} optimize {grid} {gridtmp}".split())
    shutil.move(gridtmp, grid)

    # add metadata
    metadata = source / "metadata.txt"
    runcard = next(iter(mg5_dir.glob("Events/run_01*/run_01*_tag_1_banner.txt")))
    entries = []
    if metadata.exists():
        for line in metadata.read_text().splitlines():
            k, v = line.split("=")
            entries += ["--entry", k, f"'{v}'"]
    subprocess.run(
        f"{pineappl} set {grid} {gridtmp}".split()
        + entries
        + f"--entry_from_file runcard {runcard}".split()
    )
    shutil.move(gridtmp, grid)

    # find out which PDF set was used to generate the predictions
    pdf = re.search(r"set lhaid (\d+)", (dest / "launch.txt").read_text())[1]

    # (re-)produce predictions
    with open(dest / "pineappl.convolute", "w") as fd:
        subprocess.run(
            f"{pineappl} convolute {grid} {pdf} --scales 9 --absolute --integrated".split(),
            stdout=fd,
        )
    with open(dest / "pineappl.orders", "w") as fd:
        subprocess.run(f"{pineappl} orders {grid} {pdf} --absolute".split(), stdout=fd)
    with open(dest / "pineappl.pdf_uncertainty", "w") as fd:
        subprocess.run(
            f"{pineappl} pdf_uncertainty --threads=1 {grid} {pdf}".split(), stdout=fd
        )

    return (dest / "pineappl.convolute").read_text().splitlines()[2:-2]


def results(dest, mg5_dir):
    madatnlo = next(iter(mg5_dir.glob("Events/run_01*/MADatNLO.HwU"))).read_text()
    table = filter(
        lambda line: re.match("^  [+-]", line) is not None, madatnlo.splitlines()
    )
    df = pd.DataFrame(np.array([[float(x) for x in l.split()] for l in table]))
    # start column from 1
    df.columns += 1
    df["result"] = df[3]
    df["error"] = df[4]
    df["sv_min"] = df[6]
    df["sv_max"] = df[7]

    return df


def annotate_versions(name, dest):
    grid = dest / f"{name}.pineappl"
    gridtmp = dest / f"{name}.pineappl.tmp"
    results_log = dest / "results.log"
    pineappl = paths.pineappl_exe()

    runcard_gitversion = pygit2.Repository(paths.root).describe(
        always_use_long_format=True,
        describe_strategy=pygit2.GIT_DESCRIBE_TAGS,
        dirty_suffix="-dirty",
        show_commit_oid_as_fallback=True,
    )
    mg5amc_revno = (
        subprocess.run("brz revno".split(), cwd=paths.mg5amc, stdout=subprocess.PIPE)
        .stdout.decode()
        .strip()
    )
    mg5amc_repo = (
        subprocess.run("brz info".split(), cwd=paths.mg5amc, stdout=subprocess.PIPE)
        .stdout.decode()
        .strip()
    )
    mg5amc_repo = re.search(r"\s*parent branch:\s*(.*)", mg5amc_repo)[1]

    entries = []
    entries += ["--entry", "runcard_gitversion", runcard_gitversion]
    entries += ["--entry", "mg5amc_revno", mg5amc_revno]
    entries += ["--entry", "mg5amc_repo", mg5amc_repo]
    entries += ["--entry", "lumi_id_types", "pdg_mc_ids"]
    subprocess.run(
        f"{pineappl} set {grid} {gridtmp}".split()
        + f"--entry_from_file results {results_log}".split()
        + entries
    )
    shutil.move(gridtmp, grid)


def postrun(name, dest):
    source = paths.runcards / name
    mg5_dir = dest / name
    grid = dest / f"{name}.pineappl"

    if os.access((source / "postrun.sh"), os.X_OK):
        shutil.copy2(source / "postrun.sh", dest)
        os.environ["GRID"] = str(grid)
        subprocess.run("./postrun.sh", cwd=dest)

    with lz4.frame.open(
        grid.with_suffix(grid.suffix + ".lz4"),
        "wb",
        compression_level=lz4.frame.COMPRESSIONLEVEL_MAX,
    ) as fd:
        fd.write(grid.read_bytes())
