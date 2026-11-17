import argparse, json, sys
from .core import Invalid,build,verify,cloud,read,summary

def main():
    p=argparse.ArgumentParser(description='Compare point samples, not real-world surface truth.')
    sub=p.add_subparsers(dest='command',required=True)
    b=sub.add_parser('compare');b.add_argument('contract');b.add_argument('reference');b.add_argument('reconstruction');b.add_argument('--out',required=True);b.add_argument('--require-pass',action='store_true')
    v=sub.add_parser('verify');v.add_argument('bundle');v.add_argument('--require-pass',action='store_true')
    i=sub.add_parser('inspect');i.add_argument('cloud');i.add_argument('--unit',required=True,choices=['m','mm'])
    a=p.parse_args()
    try:
        if a.command=='inspect':
            print(json.dumps({'sample_summary':summary(cloud(read(a.cloud),a.unit)),'quality_evaluated':False},sort_keys=True));return 0
        report=build(a.contract,a.reference,a.reconstruction,a.out) if a.command=='compare' else verify(a.bundle)
        print(json.dumps({'consistent':True,'status':report['status'],'precision':report['precision'],'recall':report['recall']},sort_keys=True))
        return 3 if a.require_pass and report['status']!='pass' else 0
    except Invalid as e:
        print('Invalid: '+str(e),file=sys.stderr);return 2
    except (OSError,ValueError,OverflowError):
        print('Invalid input, inconsistent bundle, or unavailable output. No path or input content is echoed.',file=sys.stderr);return 2
if __name__=='__main__':sys.exit(main())
