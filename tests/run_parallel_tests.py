
import sys
import subprocess
from pathlib import Path
from multiprocessing import cpu_count
from concurrent.futures import ProcessPoolExecutor, as_completed

def run_test_file(path, render=False):
    """Run a single test file using unittest"""

    command = [sys.executable, "-m", "unittest", str(path)]
    if render:
        command = ['RENDER_RMD=TRUE', 'RENDER_QMD=TRUE'] + command

    try:
        result = subprocess.run(
            command, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        return (path.name, result.returncode, result.stdout)
    except Exception as e:
        return (path.name, 1, f"ERROR: {e}")


def main(test_paths, max_workers=None, verbose=False, **kwargs):
    # Determine which test files to run
    if not test_paths:
        test_files = list(Path.cwd().glob("test_*.py"))
    else:
        test_files = [Path(p) for p in test_paths]

    if not test_files:
        print(f"No test files found to run.")
        sys.exit(1)

    n_cores = cpu_count() if max_workers is None else max_workers
    total = len(test_files)
    print(f"Running {total} test file(s) on {n_cores} cores...")

    passed = 0
    failed = 0
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        futures = {executor.submit(run_test_file, f, **kwargs): f for f in test_files}
        for future in as_completed(futures):
            name, code, output = future.result()
            if code == 0:
                print(f"\t✅ {name}")
                passed += 1
            else:
                print(f"\t❌ {name}")
                if verbose:
                    print(output)
                failed += 1

    # Summary
    if failed > 0:
        print(f'❌ {failed} tests failed.')
    print(f"Summary: {passed} of {total} tests passed.")

    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run unittest files in parallel.")
    parser.add_argument(
        '--render', action='store_true', default=False,
        help='Pass RENDER_RMD and RENDER_QMD arguments to the test files'
    )
    parser.add_argument(
        '--jobs', type=int,
        help='Number of parallel processes (default: # of CPU cores)'
    )
    parser.add_argument(
        '--verbose', action='store_true', default=False,
        help='Print combined output for any failing tests'
    )
    parser.add_argument(
        'tests', nargs='*',
        help='List of test files to run (default: all test_*.py in cwd)'
    )
    args = parser.parse_args()
    main(
        test_paths=args.tests,
        max_workers=args.jobs,
        render=args.render,
        verbose=args.verbose
    )
