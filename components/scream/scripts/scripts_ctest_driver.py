from utils import run_cmd, expect, check_minimum_python_version
from machines_specs import is_machine_supported

check_minimum_python_version(3, 4)

import os, shutil
from pathlib import Path

###############################################################################
class ScriptsCtestDriver(object):
###############################################################################

    ###########################################################################
    def __init__(self, submit=False, machine=None, root_dir=None, work_dir=None, dry_run=False):
    ###########################################################################

        self._submit   = submit
        self._machine  = machine
        self._root_dir = root_dir
        self._work_dir = work_dir
        self._dry_run  = dry_run

        # Probe machine if none was specified
        if self._machine is None:
            # We could potentially integrate more with CIME here to do actual
            # nodename probing.
            if "CIME_MACHINE" in os.environ and is_machine_supported(os.environ["CIME_MACHINE"]):
                self._machine = os.environ["CIME_MACHINE"]
            else:
                expect(False,
                       "scripts-ctest-driver requires either the machine arg or CIME_MACHINE in env")

        # Compute root dir (where repo is) and work dir (where build/test will happen)
        if not self._root_dir:
            self._root_dir = Path(__file__).resolve().parent.parent
        else:
            self._root_dir = Path(self._root_dir).resolve()
            expect(self._root_dir.is_dir() and self._root_dir.parts()[-2:] == ('scream', 'components'),
                   "Bad root-dir '{}', should be: $scream_repo/components/scream".format(self._root_dir))

        if self._work_dir is None:
            self._work_dir = self._root_dir.absolute().joinpath("ctest-build-scripts")
        else:
            self._work_dir = Path(self._work_dir).absolute()

        if self._work_dir.exists():
            expect(self._work_dir.is_dir(), "Work dir {} exists but is not a directory".format(self._work_dir))
            shutil.rmtree(str(self._work_dir))

        self._work_dir.mkdir(parents=True)

    ###############################################################################
    def generate_ctest_config(self, extra_configs):
    ###############################################################################
        result = "CTEST_PARALLEL_LEVEL=1 CIME_MACHINE={} ".format(self._machine)

        result += "ctest -j1 -V --output-on-failure "

        if not self._submit:
            result += "-DNO_SUBMIT=True "

        for key, value in extra_configs:
            result += "-D{}={} ".format(key, value)

        result += "-DSCREAM_ROOT={} ".format(self._root_dir)
        result += "-DBUILD_WORK_DIR={} ".format(self._work_dir)
        result += '-S {}/cmake/ctest_script_scripts_tests.cmake '.format(self._root_dir)

        return result

    ###############################################################################
    def scripts_ctest_driver(self):
    ###############################################################################
        success = True

        ctest_config = self.generate_ctest_config([])
        success = run_cmd(ctest_config, from_dir=self._work_dir, arg_stdout=None, arg_stderr=None, verbose=True, dry_run=self._dry_run)[0] == 0

        return success
