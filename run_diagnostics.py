#!/usr/bin/env python3
"""Run all bundled mathematical diagnostics in separate disposable copies.

The original research files and saved results are never used as execution
outputs. By default only a summary is printed; --report saves a detailed JSON
report to the explicitly selected path. Passing diagnostics is not a proof.
"""
from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent
RESEARCH = ROOT / 'research'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def differences(old, new, key='$'):
    if type(old) is not type(new):
        return [{'key': key, 'saved': old, 'rerun': new}]
    if isinstance(old, dict):
        result = []
        for name in sorted(set(old) | set(new)):
            result += differences(old.get(name), new.get(name), key+'.'+name)
        return result
    if isinstance(old, list):
        if len(old) != len(new):
            return [{'key': key+'.length', 'saved': len(old), 'rerun': len(new)}]
        return [change for i, (a, b) in enumerate(zip(old, new))
                for change in differences(a, b, key+'['+str(i)+']')]
    return [] if old == new else [{'key': key, 'saved': old, 'rerun': new}]

def run_one(script):
    relative = script.relative_to(RESEARCH)
    result = {'script': str(relative), 'script_sha256': digest(script)}
    with tempfile.TemporaryDirectory(prefix='sampling-diagnostic-') as directory:
        work = Path(directory) / 'research'
        shutil.copytree(RESEARCH, work)
        for helper in ('INVENTORY.json', 'publication_sources.py', 'run_frozen_diagnostic.py'):
            shutil.copyfile(ROOT/helper, Path(directory)/helper)
        shutil.copytree(ROOT/'verification', Path(directory)/'verification')
        before = {str(p.relative_to(work)): digest(p) for p in work.rglob('*.json')}
        saved_json = {str(p.relative_to(work)): json.loads(p.read_text()) for p in work.rglob('*.json')}
        env = dict(os.environ)
        env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
        start = time.monotonic()
        try:
            is_new = str(script.relative_to(ROOT)) in {e['path'] for e in json.loads((ROOT/'verification/V13-FROZEN-SELECTION.json').read_text())['files']}
            command = [sys.executable, str(Path(directory)/'run_frozen_diagnostic.py'), str(work/relative)] if is_new else [sys.executable, str(work/relative)]
            process = subprocess.run(command, cwd=work,
                env=env, capture_output=True, text=True, timeout=300)
            result.update(exit_code=process.returncode, passed=process.returncode == 0,
                status='passed' if process.returncode == 0 else 'failed',
                stdout=process.stdout, stderr=process.stderr)
            result['changed_result_files'] = [
                {'path': str(p.relative_to(work)), 'saved_sha256': before.get(str(p.relative_to(work))),
                 'rerun_sha256': digest(p),
                 'json_differences': differences(saved_json.get(str(p.relative_to(work))), json.loads(p.read_text()))}
                for p in sorted(work.rglob('*.json'))
                if before.get(str(p.relative_to(work))) != digest(p)]
        except subprocess.TimeoutExpired:
            result.update(passed=False, status='timeout_after_300_seconds')
        except OSError as error:
            result.update(passed=False, status='execution_error', error=str(error))
        result['elapsed_seconds'] = round(time.monotonic() - start, 3)
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, help='Optional detailed JSON report output path.')
    parser.add_argument('--new-only', action='store_true', help='Run only the V13 frozen additions.')
    args = parser.parse_args()
    if args.report:
        report = args.report.resolve()
        if report == ROOT or ROOT in report.parents:
            parser.error('--report must resolve outside the immutable archive directory')
    extra=json.loads((ROOT/'verification/V13-FROZEN-SELECTION.json').read_text())['diagnostic_entrypoints']
    scripts = sorted(set(RESEARCH.rglob('check_*.py')) | {ROOT/p for p in extra})
    if args.new_only:
        selected={e['path'] for e in json.loads((ROOT/'verification/V13-FROZEN-SELECTION.json').read_text())['files']}
        scripts=[p for p in scripts if str(p.relative_to(ROOT)) in selected]
    before = {str(p.relative_to(ROOT)): digest(p) for p in RESEARCH.rglob('*') if p.is_file()}
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(run_one,s) for s in scripts]
        results=[]
        for future in concurrent.futures.as_completed(futures):
            result=future.result();results.append(result)
            print(('PASS' if result['passed'] else 'FAIL')+' '+result['script'],flush=True)
        results.sort(key=lambda r:r['script'])
    for result in results:
        print(('PASS' if result['passed'] else 'FAIL') + ' ' + result['script'], flush=True)
    after = {str(p.relative_to(ROOT)): digest(p) for p in RESEARCH.rglob('*') if p.is_file()}
    summary = {
        'run_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'python': sys.version.split()[0],
        'packages': {name: importlib.metadata.version(name) for name in ('numpy', 'scipy', 'mpmath', 'sympy')},
        'script_count': len(scripts), 'passed_count': sum(r['passed'] for r in results),
        'failed_count': sum(not r['passed'] for r in results), 'not_run_count': len(scripts)-len(results),
        'all_pass': bool(scripts) and all(r['passed'] for r in results),
        'archived_research_bytes_unchanged': before == after,
        'scope': 'Diagnostic checks only. Every script ran from its own disposable copy of the complete research payload. Saved results and source bytes were not overwritten. Finite algebra and numerical fixtures do not replace analytic proofs of source or outer theorems, execute the full nonlinear sampler, or establish an efficiency improvement.',
        'results': results}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k != 'results'}, indent=2))
    return 0 if summary['all_pass'] and summary['archived_research_bytes_unchanged'] else 1

if __name__ == '__main__':
    raise SystemExit(main())
