"""Run with python3 -m unittest discover -s tests -v; no HPC tools required."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
STUB = '''#!/usr/bin/env python3
import json, os, pathlib, sys
with open(os.environ["OTB_TEST_LOG"], "a") as log:
  log.write(json.dumps({"tool": pathlib.Path(sys.argv[0]).name,
    "args": sys.argv[1:], "cwd": os.getcwd(),
    "cache": os.environ.get("NXF_SINGULARITY_CACHEDIR"),
    "library": os.environ.get("NXF_SINGULARITY_LIBRARYDIR")}) + "\\n")
'''


class ShellRegressions(unittest.TestCase):
  def setUp(self):
    self.temp = tempfile.TemporaryDirectory(prefix="otb regression ")
    self.addCleanup(self.temp.cleanup)
    self.work = Path(self.temp.name).resolve()
    for name in ("otb.sh", "run.nf"):
      shutil.copy2(ROOT / name, self.work / name)
    for name in ("scr", "config"):
      shutil.copytree(ROOT / name, self.work / name)
    self.bin = self.work / "bin"
    self.bin.mkdir()
    for name in ("nextflow", "singularity", "module"):
      self.stub(self.bin / name)
    self.log = self.work / "commands.jsonl"
    self.env = dict(os.environ)
    for name in ("NXF_SINGULARITY_CACHEDIR", "NXF_SINGULARITY_LIBRARYDIR",
                 "BUSCO", "LINEAGE", "BUSCOPATH", "POLISHTYPE", "YAHS"):
      self.env.pop(name, None)
    self.env.update(PATH=str(self.bin) + os.pathsep + self.env["PATH"],
                    OTB_TEST_LOG=str(self.log))
    (self.work / "RawData").mkdir()
    for name in ("one.bam", "two.bam"):
      (self.work / "RawData" / name).touch()

  def stub(self, path):
    path.write_text(STUB)
    path.chmod(0o755)

  def calls(self):
    return [json.loads(line) for line in self.log.read_text().splitlines()]

  def run_wrapper(self):
    result = subprocess.run(
      ["bash", "otb.sh", "--lite", "--mode", "default", "--reads",
       "RawData/*.bam", "--name", "regression", "--supress", "--check"],
      cwd=self.work, env=self.env, capture_output=True, text=True)
    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
    self.assertNotIn("NXF_SINGULARITY_CACHEDIR not set", result.stderr)
    self.assertNotIn("Nextflow Singularity cache directory set:", result.stderr)
    calls = self.calls()
    run = next(call for call in calls
               if call["tool"] == "nextflow" and call["args"][0] == "run")
    self.assertIn("--readin=RawData/*.bam", run["args"])
    container_calls = [call for call in calls if call["tool"] == "singularity"
                       and call["args"][0] in ("pull", "exec")]
    self.assertTrue(any(call["args"][0] == "pull" for call in container_calls))
    self.assertTrue(any(call["args"][0] == "exec" for call in container_calls))
    for call in container_calls:
      self.assertEqual(call["cwd"], run["cache"])
      self.assertEqual(call["cache"], run["cache"])
    return result, run

  def test_exported_cache_and_library_with_spaces(self):
    self.env["NXF_SINGULARITY_CACHEDIR"] = str(self.work / "shared cache")
    self.env["NXF_SINGULARITY_LIBRARYDIR"] = str(self.work / "shared library")
    result, run = self.run_wrapper()
    self.assertEqual(run["cache"], self.env["NXF_SINGULARITY_CACHEDIR"])
    self.assertEqual(run["library"], self.env["NXF_SINGULARITY_LIBRARYDIR"])
    self.assertNotIn("NXF_SINGULARITY_LIBRARYDIR not set", result.stderr)
    self.assertTrue(Path(run["cache"]).is_dir())

  def test_unset_cache_uses_exported_fallback(self):
    _, run = self.run_wrapper()
    self.assertEqual(run["cache"], str(self.work / "work" / "singularity"))
    self.assertTrue(Path(run["cache"]).is_dir())

  def test_library_preserved_when_cache_unset(self):
    self.env["NXF_SINGULARITY_LIBRARYDIR"] = str(self.work / "existing library")
    _, run = self.run_wrapper()
    self.assertEqual(run["library"], self.env["NXF_SINGULARITY_LIBRARYDIR"])
    self.assertEqual(run["cache"], str(self.work / "work" / "singularity"))

  def test_templates_preserve_read_patterns(self):
    self.stub(self.work / "otb.sh")
    for name in ("otb.template.slurm", "otb.lite.template.slurm"):
      for busco in ("", "--busco"):
        for reads in ("RawData/*.bam", "RawData/one.bam", "Raw Data/*.bam"):
          with self.subTest(template=name, busco=busco, reads=reads):
            script = (ROOT / name).read_text()
            script = script.replace("CCS='RawData/*.fastq'", "CCS='" + reads + "'")
            script = script.replace('Busco="--busco"', 'Busco="' + busco + '"')
            script = script.replace('Busco="" #Busco will not be run',
                                    'Busco="' + busco + '" #Busco will not be run')
            self.log.unlink(missing_ok=True)
            subprocess.run(["bash", "-c", script], cwd=self.work, env=self.env,
                           capture_output=True, text=True, check=True)
            call = next(call for call in self.calls() if call["tool"] == "otb.sh")
            args = call["args"]
            self.assertEqual(args[args.index("-in") + 1], reads)
            self.assertEqual("--busco" in args, bool(busco))
            self.assertEqual("--lite" in args, "lite" in name)

  def test_helpers_use_cache_and_explicit_location_with_spaces(self):
    cache = self.work / "shared cache"
    location = self.work / "explicit location"
    cache.mkdir()
    location.mkdir()
    for name in ("prefetch_containers.sh", "force_prefetch_containers.sh",
                 "check_containers.sh"):
      for use_cache in (False, True):
        with self.subTest(helper=name, use_cache=use_cache):
          self.log.unlink(missing_ok=True)
          if use_cache:
            self.env["NXF_SINGULARITY_CACHEDIR"] = str(cache)
          else:
            self.env.pop("NXF_SINGULARITY_CACHEDIR", None)
          subprocess.run(["bash", "scr/" + name, "-l", str(location)],
                         cwd=self.work, env=self.env, capture_output=True,
                         text=True, check=True)
          for call in self.calls():
            self.assertEqual(call["cwd"], str(cache if use_cache else location))

  def test_helpers_reject_unavailable_cache(self):
    for name in ("prefetch_containers.sh", "force_prefetch_containers.sh",
                 "check_containers.sh"):
      with self.subTest(helper=name):
        self.log.unlink(missing_ok=True)
        self.env["NXF_SINGULARITY_CACHEDIR"] = str(self.work / "missing cache")
        result = subprocess.run(["bash", "scr/" + name], cwd=self.work,
                                env=self.env, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unable to enter Singularity cache directory", result.stderr)
        self.assertFalse(self.log.exists())


if __name__ == "__main__":
  unittest.main()
