"""Workshop-only parameter sensitivity and target matching on saved energies."""
import csv
import hashlib
import json
import math
from pathlib import Path
import webbrowser
from source_profile import read_qfep,reconstruct_source
from kemp_profile import write_outputs


def published_reference(case):
    values={'Des27_7_IN':(13.2,.6),'Des27_7_OUT':(15.5,1.5),'Des27_WT_IN':(16.6,2.1),'Des27_WT_OUT':(17.1,2.4)}
    if case not in values:return None
    mean,sd=values[case]
    return {'case':case,'barrier_kcal_mol':mean,'replica_sd_kcal_mol':sd,'doi':'10.1038/s41586-025-09136-2','location':'Supplementary Table3','interpretation':'Reported mean/SD, not a full paper profile.'}


def scan(raw,couplings,offsets,target=13.2,bins=18,min_points=20,skip=0,refine=True):
    if not math.isfinite(target) or target<=0:raise ValueError('Positive finite workshop target required.')
    hs=sorted(set(couplings));alphas=sorted(set(offsets))
    if not hs or not alphas or any(not math.isfinite(h) or h<=0 for h in hs) or any(not math.isfinite(a) for a in alphas):raise ValueError('Use finite positive couplings and finite offsets.')
    if len(hs)*len(alphas)>2500:raise ValueError('Limit a workshop grid to2500 combinations.')
    archived=read_qfep();cache={}
    def evaluate(h,a,kind='grid'):
        key=(float(h),float(a))
        if key not in cache:
            trial={**archived,'coupling_kcal_mol':float(h),'native_alpha_kcal_mol':float(a)}
            analysis,summary=reconstruct_source(raw,bins,min_points,skip,trial)
            is_archived=abs(h-archived['coupling_kcal_mol'])<1e-10 and abs(a-archived['native_alpha_kcal_mol'])<1e-10
            summary['native_qfep']=archived
            summary.update({'workshop_only':True,'candidate_is_archived':is_archived,
                'parameter_origin':'archived comparison' if is_archived else 'workshop candidate, not archived parameters',
                'plot_title':'Kemp EVB — WORKSHOP PARAMETER TRIAL; NOT A PREDICTION',
                'reference':'reactant' if analysis['thermochemistry']['reactant'] else 'minimum'})
            summary['warnings'][0]='Workshop sensitivity/target-matching exercise; candidate parameters are not independent physical calibration.'
            barrier=summary['conditional_barrier_estimate_kcal_mol'];dg=summary['conditional_reaction_free_energy_kcal_mol']
            curve=analysis['profile']['state_free_energies_kcal_mol']['evb_ground']
            reactant=analysis['thermochemistry']['reactant'];zero=reactant['free_energy_kcal_mol'] if reactant else min((v for v in curve if v is not None),default=0.)
            profile_peak=None
            transitions=analysis['thermochemistry']['transition_states']
            if barrier is not None and transitions:
                peak=max(transitions,key=lambda item:item['free_energy_kcal_mol']);i=int(peak['index'])
                centres=analysis['profile']['gap_centers_kcal_mol'];adjusted=float(centres[i])
                raw_gap=adjusted+summary['total_product_offset_kcal_mol']
                width=float(centres[1]-centres[0])
                profile_peak={'bin_index':i,'raw_gap_kcal_mol':raw_gap,'offset_adjusted_gap_kcal_mol':adjusted,
                    'free_energy_relative_reactant_kcal_mol':float(curve[i]-zero),'bin_width_kcal_mol':width,
                    'raw_bin_lower_kcal_mol':raw_gap-width/2,'raw_bin_upper_kcal_mol':raw_gap+width/2,
                    'raw_peak_bin_contains_zero':bool(raw_gap-width/2<=0<=raw_gap+width/2),
                    'offset_adjusted_peak_bin_contains_zero':bool(adjusted-width/2<=0<=adjusted+width/2),
                    'interpretation':'Binned free-energy maximum (TS estimate); a molecular transition state has not been validated.'}
            record={'id':f'candidate_{len(cache)}','kind':kind,'coupling_kcal_mol':float(h),'native_alpha_kcal_mol':float(a),
                'total_product_offset_kcal_mol':summary['total_product_offset_kcal_mol'],'barrier_kcal_mol':barrier,'reaction_dg_kcal_mol':dg,
                'target_error_kcal_mol':None if barrier is None else abs(barrier-target),'internal_missing_bins':summary['internal_missing_bins'],
                'candidate_is_archived':is_archived,'reference':summary['reference'],'profile_peak':profile_peak,'x':[float(x)+summary['total_product_offset_kcal_mol'] for x in analysis['profile']['gap_centers_kcal_mol']],
                'y':[None if v is None else float(v-zero) for v in curve]}
            cache[key]=(record,analysis,summary)
        return cache[key]
    baseline=evaluate(archived['coupling_kcal_mol'],archived['native_alpha_kcal_mol'],'archived')[0]
    grid=[]
    for a in alphas:
        for h in hs:grid.append(evaluate(h,a)[0])
    matches=[]
    if refine:
        for a in alphas:
            found=[]
            for left,right in zip(hs,hs[1:]):
                l=evaluate(left,a)[0];r=evaluate(right,a)[0]
                if l['barrier_kcal_mol'] is None or r['barrier_kcal_mol'] is None:continue
                fl=l['barrier_kcal_mol']-target;fr=r['barrier_kcal_mol']-target
                if fl*fr>0:continue
                lo,hi=left,right
                for _ in range(45):
                    mid=(lo+hi)/2;row=evaluate(mid,a,'refinement')[0];value=row['barrier_kcal_mol']
                    if value is None:break
                    fm=value-target
                    if abs(fm)<1e-5:
                        if not any(abs(mid-h)<1e-5 for h in found):
                            row['kind']='target-matched';matches.append(row);found.append(mid)
                        break
                    if fl*fm<=0:hi=mid
                    else:lo=mid;fl=fm
    valid=[row for row in grid if row['barrier_kcal_mol'] is not None]
    if not valid:raise ValueError('No grid candidate has a resolved barrier. Retain gaps; no peak will be inserted. Inspect bin counts/sampling.')
    if matches:
        best=min(matches,key=lambda row:(abs(row['native_alpha_kcal_mol']-archived['native_alpha_kcal_mol']),abs(row['coupling_kcal_mol']-archived['coupling_kcal_mol']),row['target_error_kcal_mol']))
    else:best=min(valid,key=lambda row:row['target_error_kcal_mol'])
    record,analysis,summary=evaluate(best['coupling_kcal_mol'],best['native_alpha_kcal_mol'])
    summary.update({'profile_peak':best['profile_peak'],'workshop_target_kcal_mol':target,'target_matching_error_kcal_mol':best['target_error_kcal_mol'],
        'selection_policy':'Among target matches, prefer the scanned offset closest to the archived offset, then coupling closest to the archived coupling; otherwise choose closest grid barrier.',
        'fit_type':'coupling target-match at a fixed scanned offset' if best['kind']=='target-matched' else 'closest workshop grid candidate'})
    if best['kind']=='target-matched':
        summary['target_barrier_kcal_mol']=target
        summary['plot_title']='Kemp EVB — WORKSHOP TARGET MATCH; NOT A PREDICTION'
    displayed=grid+matches
    if baseline not in displayed:displayed.append(baseline)
    return {'workshop_only':True,'calibrated':False,'convergence_established':False,'physical_barrier_available':False,
        'case':raw.get('metadata',{}).get('case'),'published_reference':published_reference(raw.get('metadata',{}).get('case')),
        'target_barrier_kcal_mol':target,'bins':bins,'minimum_per_window_bin':min_points,'extra_skip_samples':skip,
        'archived_parameters':archived,'baseline':baseline,'grid_candidates':grid,'target_matches':matches,'displayed_candidates':displayed,
        'best':best,'numerical_evaluations':len(cache),'grid_couplings':hs,'grid_native_offsets':alphas,
        'warnings':['One barrier does not uniquely determine coupling and offset.','Target matches are workshop demonstrations, not independently predicted enzyme barriers.',
            'No new sampling, interpolation, synthetic peak or vertical energy rescaling is used.','Same short correlated trajectory and changed-model validation limits apply.']},analysis,summary


def write_scan(result,analysis,summary,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/'scan.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    fields=['kind','coupling_kcal_mol','native_alpha_kcal_mol','total_product_offset_kcal_mol','barrier_kcal_mol','reaction_dg_kcal_mol','target_error_kcal_mol','candidate_is_archived','profile_peak_raw_gap_kcal_mol','profile_peak_offset_adjusted_gap_kcal_mol','profile_peak_bin_width_kcal_mol','raw_peak_bin_contains_zero']
    with (out/'combinations.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,extrasaction='ignore');writer.writeheader()
        for row in sorted(result['displayed_candidates'],key=lambda r:float('inf') if r['target_error_kcal_mol'] is None else r['target_error_kcal_mol']):
            peak=row['profile_peak'] or {}
            writer.writerow({**row,'profile_peak_raw_gap_kcal_mol':peak.get('raw_gap_kcal_mol'),
                'profile_peak_offset_adjusted_gap_kcal_mol':peak.get('offset_adjusted_gap_kcal_mol'),
                'profile_peak_bin_width_kcal_mol':peak.get('bin_width_kcal_mol'),'raw_peak_bin_contains_zero':peak.get('raw_peak_bin_contains_zero')})
    write_outputs(analysis,summary,out/'selected_example')
    payload=json.dumps(result,allow_nan=False).replace('<','\\u003c')
    page=Path(__file__).with_name('workshop_viewer.html').read_text().replace('__SCAN_DATA__',payload)
    (out/'scan.html').write_text(page)


def run(args,out,seed):
    path=Path(out)/f'replica_{seed}_raw.json';raw=json.loads(path.read_text())
    target=13.2 if args.target_barrier is None else args.target_barrier
    result,analysis,summary=scan(raw,args.couplings,args.offsets,target,args.bins,args.min_points,args.skip_samples,not args.no_refine)
    result['raw_sha256']=hashlib.sha256(path.read_bytes()).hexdigest();summary['raw_sha256']=result['raw_sha256']
    folder=Path(out)/f'workshop_scan_target{target:g}_bins{args.bins}_min{args.min_points}'
    write_scan(result,analysis,summary,folder)
    print(f'Workshop scan: {folder/"scan.html"}')
    print(f'{len(result["grid_candidates"])} grid combinations; {len(result["target_matches"])} target-matched examples. No new dynamics.')
    b=result['best'];print(f'Selected workshop example: H12={b["coupling_kcal_mol"]:.6f}; alpha={b["native_alpha_kcal_mol"]:.6f}; barrier={b["barrier_kcal_mol"]:.6f}; reaction DG={b["reaction_dg_kcal_mol"]:.6f} kcal/mol. Target-matched, not independently predicted.')
    if args.open:webbrowser.open((folder/'scan.html').resolve().as_uri())
