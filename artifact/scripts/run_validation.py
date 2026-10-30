"""Run documented checks, retaining one log per command and machine-readable status."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np

JOBS = [
    ('bridge_realizability', ['bridge_realizability_certificate.py']),
    ('target_graph', ['target_transition_graph.py']),
    ('binary_grammar', ['binary_target_grammar_certificate.py']),
    ('a5', ['a5_criterion.py']),
    ('bridge_tests', ['-m','unittest','test_bridge_sweeps']),
    ('derived_small', ['derived_length_test.py','small']),
    ('geometric', ['geometric_lemma.py']),
    ('gl2z4', ['gl2z4.py']),
    ('noninvolution', ['noninvolution_leaves.py']),
    ('bridge_witnesses', ['bridge_verify_counterexample.py','--case','all']),
    ('generates', ['generates_G.py']),
    ('gl3f2', ['gl3f2.py']),
    ('pgl7_sample', ['mixing_level_pgl.py','7']),
    ('pgl7_bottleneck', ['bottleneck_pgl7.py']),
    ('strengthening', ['test_strengthening.py']),
    ('quantitative_families', ['test_quantitative_families.py']),
    ('octahedral_pair', ['octahedral_pair_check.py']),
    ('octahedral_cap6', ['octahedral_cap6_check.py']),
    ('hamilton_reduction', ['s4_hamilton_reduction_check.py']),
    ('s4_internal_gap', ['s4_internal_gap_certify.py']),
    ('ordered_groups', ['verification_groups.py']),
    ('rank_one_witnesses', ['pgl2_classification.py','5','7','9','11','13','17','--tries','100']),
    ('psl33', ['psl33_maximals_check.py']),
    ('transvections', ['transvection_levels.py']),
    ('allpairs_evidence', ['../external/p2a/s4_allpairs_doublecoset_audit.py','11','13','19','29','37']),
    ('p0_external', ['../external/p0/20_level3_exact_certificate.py']),
]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--jobs',type=int,default=2)
    args=parser.parse_args()
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=True)
    scripts=Path(__file__).resolve().parent
    env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
             MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', PYTHONHASHSEED='0',
             PYTHONUNBUFFERED='1')
    def run(job):
        name,argv=job; begin=time.monotonic()
        with (output/(name+'.log')).open('w') as log:
            r=subprocess.run([sys.executable]+argv,cwd=scripts,env=env,
                             stdout=log,stderr=subprocess.STDOUT)
        return dict(name=name,arguments=argv,returncode=r.returncode,
                    seconds=round(time.monotonic()-begin,3),
                    log_sha256=hashlib.sha256((output/(name+'.log')).read_bytes()).hexdigest())
    rows=[]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for future in as_completed([pool.submit(run,j) for j in JOBS]):
            row=future.result(); rows.append(row)
            print(row['name'],row['returncode'],row['seconds'],flush=True)
            report=dict(python=platform.python_version(),numpy=np.__version__,
                        results=sorted(rows,key=lambda r:r['name']))
            (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return int(any(r['returncode'] for r in rows))


if __name__=='__main__':
    raise SystemExit(main())
