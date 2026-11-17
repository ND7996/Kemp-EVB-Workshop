"""Offline, parameter-conditional Kemp profile reconstruction (no new dynamics)."""
import csv
import html
import json
import math
from pathlib import Path
import webbrowser
from evbopenmm.models import FEPSampleRecord,FEPWindowRecord
from evbopenmm.analysis.mapping_fep import analyze_mapping_fep_umbrella_records
from evbopenmm.analysis.profile import analyze_evb_profile_records


def reconstruct(raw,coupling,extra_offset,bins,min_points):
    if raw.get('schema')!='kemp_full_mapping_records_v1':raise ValueError('Expected full-enzyme mapping records, not diagnostics or the obsolete trial schema.')
    if not math.isfinite(coupling) or coupling<=0 or not math.isfinite(extra_offset):raise ValueError('Coupling must be finite and positive; offset must be finite.')
    if bins<2 or min_points<1:raise ValueError('Use at least two bins and one point per bin.')
    correction=float(raw['metadata']['morse_P_minus_R_correction_kcal_mol'])
    windows=[]
    for window in raw['windows']:
        lam=float(window['lambda_R']);samples=[]
        for row in window['samples']:
            r,p=float(row['E_R_kj_mol']),float(row['E_P_kj_mol'])
            if not all(math.isfinite(v) for v in (r,p)):raise ValueError('Nonfinite raw energy.')
            samples.append(FEPSampleRecord(int(row['step']),lam*r+(1-lam)*p,r-p,(r,p)))
        if not samples:raise ValueError('Each selected mapping window needs samples; complete or resume sampling first.')
        windows.append(FEPWindowRecord(int(window['index']),(lam,1-lam),tuple(samples)))
    temperature=float(raw['config']['temperature_kelvin'])
    base=analyze_mapping_fep_umbrella_records(windows,state_names=('R','P'),kT_kj_mol=.00831446261815324*temperature,
        skip_samples=0,n_gap_bins=bins,min_points_per_bin=min_points,input_energy_unit='kJ/mol')
    analysis=analyze_evb_profile_records(windows,state_names=('R','P'),base_analysis=base,
        diabatic_offsets_kcal_mol=(0.,correction+extra_offset),coupling_matrix_kcal_mol=((0.,coupling),(coupling,0.)),
        skip_samples=0,n_gap_bins=bins,min_points_per_bin=min_points,input_energy_unit='kJ/mol')
    profile=analysis['profile'];ground=profile['state_free_energies_kcal_mol']['evb_ground']
    valid=[i for i,v in enumerate(ground) if v is not None and math.isfinite(v)]
    missing=[i for i in range(min(valid),max(valid)+1) if i not in valid] if valid else []
    summary={'calibrated':False,'convergence_established':False,'physical_barrier_available':False,
        'coupling_kcal_mol':coupling,'extra_product_offset_kcal_mol':extra_offset,'morse_convention_correction_kcal_mol':correction,
        'total_product_offset_kcal_mol':correction+extra_offset,'resolved_bins':len(valid),'internal_missing_bins':missing,
        'conditional_barrier_estimate_kcal_mol':analysis['dg_dagger_kcal_mol'],
        'conditional_reaction_free_energy_kcal_mol':analysis['dg0_kcal_mol'],
        'profile_shape':analysis['profile_shape'],'coverage':analysis['gap_coverage'],
        'sampling_complete':raw['cursor']['phase']=='complete',
        'warnings':list(analysis['thermochemistry']['warnings'])+['Uncalibrated trial parameters; not the paper barrier or a validated physical prediction.','A single short trajectory does not establish convergence or independent-sample uncertainty.']}
    if missing:summary['warnings'].append('Internal unresolved bins: do not interpret a barrier across the disconnected profile.')
    return analysis,summary


def write_outputs(analysis,summary,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    x=list(analysis['profile']['gap_centers_kcal_mol']);y=list(analysis['profile']['state_free_energies_kcal_mol']['evb_ground'])
    finite=[v for v in y if v is not None and math.isfinite(v)];zero=min(finite) if finite else 0.
    reference=summary.get('reference','minimum')
    if reference=='reactant':zero=float(analysis['thermochemistry']['reactant']['free_energy_kcal_mol'])
    y=[None if v is None or not math.isfinite(v) else v-zero for v in y]
    with (out/'profile.csv').open('w',newline='') as handle:
        writer=csv.writer(handle);writer.writerow([summary.get('gap_csv_column','raw_energy_gap_kcal_mol'),f'trial_free_energy_relative_to_{reference}_kcal_mol','bin_count'])
        writer.writerows(zip(x,y,analysis['profile']['bin_counts']))
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    xmin,xmax=min(x),max(x);yvalues=[v for v in y if v is not None]+[0.]
    paper_refs={'Des27_7_IN':(13.2,0.6),'Des27_7_OUT':(15.5,1.5),'Des27_WT_IN':(16.6,2.1),'Des27_WT_OUT':(17.1,2.4)}
    paper=paper_refs.get(summary.get('case'))
    if paper: yvalues.extend([paper[0]-paper[1],paper[0]+paper[1]])
    span=max(max(yvalues)-min(yvalues),1.);ymin=min(yvalues)-.05*span;ymax=max(yvalues)+.08*span
    def px(v):return 90+760*(v-xmin)/max(xmax-xmin,1e-12)
    def py(v):return 420-320*(v-ymin)/(ymax-ymin)
    title='Kemp EVB — FITTED TO CHOSEN TARGET; NOT A PREDICTION' if 'target_barrier_kcal_mol' in summary else 'Kemp EVB — UNCALIBRATED TRIAL PROFILE'
    title=summary.get('plot_title',title)
    elements=['<svg xmlns="http://www.w3.org/2000/svg" width="920" height="540" viewBox="0 0 920 540">',
        '<rect width="920" height="540" fill="white"/><g font-family="sans-serif" fill="#223">',
        f'<text x="460" y="30" text-anchor="middle" font-size="19">{title}</text>',
        '<path d="M90 85V420H850" fill="none" stroke="#333"/>']
    for i in range(6):
        xv=xmin+(xmax-xmin)*i/5;yv=ymin+(ymax-ymin)*i/5
        elements.extend([f'<text x="{px(xv):.2f}" y="445" text-anchor="middle">{xv:.1f}</text>',f'<text x="80" y="{py(yv)+5:.2f}" text-anchor="end">{yv:.1f}</text>'])
    if paper:
        mu,sd=paper
        elements.extend([f'<rect x="90" y="{py(mu+sd):.2f}" width="760" height="{py(mu-sd)-py(mu+sd):.2f}" fill="#ddd" opacity="0.6"/>',f'<path d="M90 {py(mu):.2f}H850" stroke="#777" stroke-dasharray="6 4"/>',f'<text x="95" y="65" font-size="12">Blue: computed OpenMM profile. Grey: published Q6 barrier {mu:g} ± {sd:g} kcal/mol (replica SD).</text>'])
    segment=[]
    for xv,yv in zip(x,y):
        if yv is None:
            if segment:elements.append('<polyline points="'+' '.join(segment)+'" fill="none" stroke="#1769aa" stroke-width="2"/>');segment=[]
        else:
            segment.append(f'{px(xv):.2f},{py(yv):.2f}')
            elements.append(f'<circle cx="{px(xv):.2f}" cy="{py(yv):.2f}" r="3" fill="#1769aa"/>')
    if segment:elements.append('<polyline points="'+' '.join(segment)+'" fill="none" stroke="#1769aa" stroke-width="2"/>')
    if 'target_barrier_kcal_mol' in summary:
        value=summary['conditional_barrier_estimate_kcal_mol']
        elements.append(f'<path d="M90 {py(value):.2f}H850" stroke="#b33" stroke-dasharray="5 5"/><text x="500" y="{py(value)-8:.2f}" fill="#a22">Fitted forward barrier = {value:.4f} kcal/mol</text>')
    elements.extend([f'<text x="460" y="480" text-anchor="middle">{html.escape(summary.get("gap_axis_label","Energy gap (kcal/mol)"))}</text>',
        f'<text transform="translate(22 250) rotate(-90)" text-anchor="middle">Trial free energy relative to {reference} (kcal/mol)</text>',
        f'<text x="460" y="510" text-anchor="middle">H12={summary["coupling_kcal_mol"]:g}; extra offset={summary["extra_product_offset_kcal_mol"]:g} kcal/mol. Gaps are not interpolated.</text>',
        '</g></svg>'])
    svg=''.join(elements);(out/'profile.svg').write_text(svg)
    page=('<!doctype html><html lang="en"><head><meta charset="utf-8">'
          '<meta name="viewport" content="width=device-width, initial-scale=1">'
          '<title>Kemp EVB profile</title><style>body{margin:0;padding:20px;background:#f7f8f6}'
          'svg{display:block;width:100%;height:auto;max-width:1100px;margin:auto;background:white}'
          '</style></head><body>'+svg+'</body></html>')
    (out/'profile.html').write_text(page)


def run(args,out,seed):
    # All user-facing profiles use frame-wise mixing/reweighting, including
    # explicit coupling/offset trials. The old binned-curve helper remains only
    # for historical tests and explicitly labelled target-fit demonstrations.
    if args.coupling is None and args.offset!=0:
        raise ValueError('Omit --offset for the archived comparison, or provide --coupling for a labelled user-parameter trial.')
    from source_profile import run as run_source
    return run_source(args,out,seed)


def fit_target(raw,coupling,target,bins,min_points):
    """Fit only the extra P offset, holding the explicitly chosen coupling fixed.

    A barrier alone cannot identify both parameters. This is a target-matching
    exercise, not independent reference calibration or paper reproduction.
    """
    import warnings
    from evbopenmm.analysis.profile import build_diabatic_profile,analyze_evb_ground
    if not math.isfinite(target) or target<=0:raise ValueError('Target barrier must be finite and positive.')
    initial,baseline=reconstruct(raw,coupling,0.,bins,min_points)
    if baseline['internal_missing_bins'] or baseline['resolved_bins']<3:
        raise ValueError('Target fitting requires a connected resolved profile. The original 40-bin user curve is fragmented; an explicitly coarser --bins 18 can be explored without interpolation, but is not new sampling.')
    sampling=build_diabatic_profile(initial['windows'],state_names=('R','P'),base_analysis=initial,
        n_gap_bins=bins,min_points_per_bin=min_points,input_energy_unit='kJ/mol')
    correction=baseline['morse_convention_correction_kcal_mol']
    def evaluate(offset):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            result=analyze_evb_ground(sampling,diabatic_offsets_kcal_mol=(0.,correction+offset),coupling_matrix_kcal_mol=((0.,coupling),(coupling,0.)))
        return result['dg_dagger_kcal_mol']
    roots=[];previous=None
    for offset in range(-400,101,2):
        barrier=evaluate(offset)
        if barrier is not None and math.isfinite(barrier):
            residual=barrier-target
            if previous is not None and previous[1]*residual<=0:
                lo,flo=previous;hi=float(offset)
                for _ in range(60):
                    mid=.5*(lo+hi);value=evaluate(mid)
                    if value is None:break
                    fm=value-target
                    if abs(fm)<1e-6:
                        roots.append((mid,value));break
                    if flo*fm<=0:hi=mid
                    else:lo,flo=mid,fm
            previous=(float(offset),residual)
        else:previous=None
    if not roots:raise ValueError('No continuous target-matching offset found in the declared [-400,100] kcal/mol search interval. No artificial peak was added.')
    offset,_=min(roots,key=lambda item:abs(item[0]))
    analysis,summary=reconstruct(raw,coupling,offset,bins,min_points)
    achieved=summary['conditional_barrier_estimate_kcal_mol']
    if achieved is None or abs(achieved-target)>1e-4:raise ValueError('Final reconstruction did not reproduce the target; refusing to report a fitted result.')
    summary.update({'fit_type':'single-target offset fit at fixed illustrative coupling','target_barrier_kcal_mol':target,
        'target_source':'User-selected target; 13.2 kcal/mol corresponds to published Des27.7-IN computed mean, not an independent validation datum.',
        'target_residual_kcal_mol':achieved-target,'reference':'reactant','calibration_validated':False,
        'fit_offset_search_interval_kcal_mol':[-400,100],'target_matching_offsets_found':len(roots),
        'original_zero_extra_offset_barrier_at_same_binning_kcal_mol':baseline['conditional_barrier_estimate_kcal_mol']})
    summary['warnings'].extend(['FITTED TO A CHOSEN TARGET: the achieved barrier is not an independent prediction.',
        'Holding coupling fixed does not establish it physically; another coupling can give another matching offset.',
        'Coarse connected binning does not repair unresolved finer-scale sampling; bin-width and replica sensitivity remain untested.',
        'Reaction free energy and fitted state offset have no independent calibration here.'])
    return analysis,summary


def run_target(args,out,seed):
    if args.coupling is None or args.target_barrier is None:raise ValueError('Supply --coupling and --target-barrier explicitly for a labelled target fit.')
    raw=json.loads((getattr(args,'records',None) or Path(out)/f'replica_{seed}_raw.json').read_text())
    analysis,summary=fit_target(raw,args.coupling,args.target_barrier,args.bins,args.min_points)
    import hashlib
    summary['raw_sha256']=hashlib.sha256((Path(out)/f'replica_{seed}_raw.json').read_bytes()).hexdigest()
    target=Path(out)/f'target_fit_{args.target_barrier:g}_H{args.coupling:g}_bins{args.bins}_min{args.min_points}'
    write_outputs(analysis,summary,target)
    print(f'TARGET-FITTED illustration: {target / "profile.html"}')
    print(f'Fitted barrier: {summary["conditional_barrier_estimate_kcal_mol"]:.6f} kcal/mol; imposed target, not a prediction.')
    print(f'Fixed coupling={args.coupling:g}; fitted extra P offset={summary["extra_product_offset_kcal_mol"]:.6f}; conditional reaction free energy={summary["conditional_reaction_free_energy_kcal_mol"]:.6f} kcal/mol.')
    if args.open:webbrowser.open((target/'profile.html').resolve().as_uri())
