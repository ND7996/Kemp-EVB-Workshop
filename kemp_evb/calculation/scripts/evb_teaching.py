"""Reduced five-atom Kemp teaching Hamiltonian in OpenMM, not a full enzyme model.

Uses three published Q6 Morse terms and state-specific C-N harmonic terms.
Coupling and extra state offset are illustrative controls, not fitted parameters.
"""
import json
import math
import os
from pathlib import Path
import sys
import numpy as np
import openmm as mm
from openmm import app,unit as u
from kemp_inputs import ROOT,audit_case,case_paths,sections,read


def library_energy():
    if os.environ.get('EVBOPENMM_SRC'):
        p=Path(os.environ['EVBOPENMM_SRC']).expanduser()
        if (p/'src/evbopenmm').is_dir():p=p/'src'
        sys.path.insert(0,str(p))
    try:
        from evbopenmm.analysis.energetics import evb_ground_state_energy_kj_mol
    except ImportError as e:
        raise ImportError('Install your EVBOpenMM checkout with python -m pip install -e /path/to/EVBOpenMM') from e
    return evb_ground_state_energy_kj_mol


def build_core(case,coupling,offset):
    if not math.isfinite(coupling) or coupling<=0 or not math.isfinite(offset):raise ValueError('Coupling must be positive; offset must be finite.')
    data=audit_case(case);atoms={a['q_index']:a for a in data['reactive_atom_identities']}
    roles=[1,13,4,5,21] # CH, HT, OR, NR, catalytic Asp OD1
    xyz=np.array([atoms[i]['xyz_angstrom'] for i in roles])/10
    f=sections(read(case_paths(case)[0]));types={int(row[0]):list(map(float,row[1:])) for row in f['bond_types']}
    system=mm.System()
    for mass in (12.01,1.008,16.,14.007,16.):system.addParticle(mass)
    def diagonal(product):
        if product:
            d,a,r=types[4];k,rcn=types[8]
            expression=f'{d*4.184}*(exp(-2*{a*10}*(distance(p5,p2)-{r/10}))-2*exp(-{a*10}*(distance(p5,p2)-{r/10})))+0.5*{k*418.4}*(distance(p1,p4)-{rcn/10})^2'
        else:
            d,a,r=types[9];e,b,s=types[10];k,rcn=types[7]
            expression=f'{d*4.184}*(exp(-2*{a*10}*(distance(p1,p2)-{r/10}))-2*exp(-{a*10}*(distance(p1,p2)-{r/10})))+{e*4.184}*(exp(-2*{b*10}*(distance(p3,p4)-{s/10}))-2*exp(-{b*10}*(distance(p3,p4)-{s/10})))+0.5*{k*418.4}*(distance(p1,p4)-{rcn/10})^2'
        force=mm.CustomCompoundBondForce(5,expression);force.addBond(range(5),[]);return force
    lower=mm.CustomCVForce('0.5*(VR+VP+alpha-sqrt((VR-VP-alpha)^2+4*H12^2))')
    lower.addGlobalParameter('alpha',offset*4.184);lower.addGlobalParameter('H12',coupling*4.184)
    lower.addCollectiveVariable('VR',diagonal(False));lower.addCollectiveVariable('VP',diagonal(True));system.addForce(lower)
    return system,lower,xyz,data


def geometries(xyz):
    donor,_,oxygen,nitrogen,acceptor=xyz
    co=(acceptor-donor)/np.linalg.norm(acceptor-donor)
    cn=(nitrogen-donor)/np.linalg.norm(nitrogen-donor)
    opening=(oxygen-nitrogen)/np.linalg.norm(oxygen-nitrogen)
    h_start=donor+co*.108;h_end=acceptor-co*.0945
    for t in np.linspace(0,1,81):
        p=xyz.copy();p[1]=(1-t)*h_start+t*h_end
        p[3]=donor+cn*((1-t)*.132+t*.1157)
        p[2]=oxygen+opening*.14*t
        yield float(t),p


def run(case,output,coupling=10.,offset=0.):
    energy=library_energy();system,cv,xyz,data=build_core(case,coupling,offset)
    output.mkdir(parents=True,exist_ok=False)
    integ=mm.VerletIntegrator(.00025);ctx=mm.Context(system,integ,mm.Platform.getPlatformByName('Reference'))
    rows=[]
    for t,p in geometries(xyz):
        ctx.setPositions(p*u.nanometer);vr,vp=map(float,cv.getCollectiveVariableValues(ctx))
        state=ctx.getState(getEnergy=True,getForces=True)
        e=float(state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole))
        reference=energy(vr,vp,alpha_kj_mol=offset*4.184,hij_kj_mol=coupling*4.184)
        forces=np.asarray(state.getForces(asNumpy=True).value_in_unit(u.kilojoule_per_mole/u.nanometer))
        if abs(e-reference)>1e-8 or not np.isfinite(forces).all():raise RuntimeError('OpenMM and EVBOpenMM reduced-model energy mismatch.')
        gap=vp+offset*4.184-vr;p_product=.5*(1-gap/math.sqrt(gap*gap+4*(coupling*4.184)**2))
        rows.append({'path_progress':t,'VR_kcal_mol':vr/4.184,'VP_unshifted_kcal_mol':vp/4.184,
                     'Eminus_kcal_mol':e/4.184,'product_state_weight':p_product})
    summary={'status':'reduced_openmm_evb_teaching_scan_completed','case':case,'atoms':5,'path_points':len(rows),
       'coupling_kcal_mol':coupling,'extra_state_offset_kcal_mol':offset,'coupling_and_offset_are_illustrative':True,
       'included_terms':'Published three Morse and two alternative C-N harmonic terms; full two-state lower eigenvalue evaluated by OpenMM.',
       'omitted_terms':'All other bonded and nonbonded terms, protein, solvent, SCAAS, LRF, entropy and actual reaction-path optimization.',
       'geometry_path':'User-defined interpolation of proton transfer, N-O opening and C-N shortening, not a minimum-energy/free-energy path or MD.',
       'not_full_paper_reproduction':True,'barrier_kcal_mol':None,
       'validation':'OpenMM lower-surface energies checked at all 81 points against uploaded EVBOpenMM analysis.energetics; forces finite.'}
    (output/'teaching_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    # Offline interactive graphic from actual evaluated diabatic terms. Controls
    # change only illustrative alpha/H12; no molecular simulation is claimed.
    template='''<!doctype html><meta charset="utf-8"><title>Kemp EVB teaching model</title>
<style>body{font:17px sans-serif;max-width:1000px;margin:2rem auto;padding:1rem}canvas{width:100%;border:1px solid #ddd}.notice{background:#fff0cd;padding:1rem}label{display:block;margin:1rem}</style>
<h1>One EVB model, different simulation engines</h1><p class="notice">REDUCED FIVE-ATOM TEACHING MODEL. Not the paper's enzyme free-energy profile. Protein/solvent and most force-field terms are omitted. The geometric path is imposed. H12 and offset are illustrative, not fitted.</p>
<label>Coupling H12 (kcal/mol): <input id="h" type="range" min="0.1" max="40" step="0.1" value="HVALUE"><span id="hv"></span></label>
<label>Extra product-state offset (kcal/mol): <input id="a" type="range" min="-80" max="80" step="1" value="AVALUE"><span id="av"></span></label>
<canvas id="c" width="1000" height="480"></canvas><p>Blue: reactant diabat. Orange: shifted product diabat. Green: lower EVB eigenvalue. Common plot energy zero is the lowest currently displayed value. X: imposed geometric progress, not time or sampled free-energy coordinate.</p>
<p>Native Morse wells have minima -De. Changing the energy-zero convention requires explicit offset bookkeeping. These controls demonstrate mixing; they do not calibrate this enzyme.</p>
<script>const rows=DATA;const c=document.getElementById('c'),ctx=c.getContext('2d');function draw(){let H=+document.getElementById('h').value,A=+document.getElementById('a').value;document.getElementById('hv').textContent=H;document.getElementById('av').textContent=A;let r=rows.map(x=>x.VR_kcal_mol),p=rows.map(x=>x.VP_unshifted_kcal_mol+A),e=r.map((v,i)=>(v+p[i]-Math.sqrt((v-p[i])**2+4*H*H))/2);let lo=Math.min(...r,...p,...e),hi=Math.max(...r,...p,...e);ctx.clearRect(0,0,1000,480);ctx.strokeStyle='#444';ctx.beginPath();ctx.moveTo(70,30);ctx.lineTo(70,410);ctx.lineTo(965,410);ctx.stroke();[r,p,e].forEach((v,k)=>{ctx.strokeStyle=['#1769aa','#df7926','#34884d'][k];ctx.lineWidth=3;ctx.beginPath();v.forEach((y,i)=>{let x=70+i/(v.length-1)*895,z=410-(y-lo)/(hi-lo)*360;i?ctx.lineTo(x,z):ctx.moveTo(x,z)});ctx.stroke()});ctx.fillStyle='#333';ctx.font='16px sans-serif';ctx.fillText('0',65,435);ctx.fillText('Imposed geometric progress',380,455);ctx.fillText('1',950,435);ctx.fillText('Potential energy, not free energy (kcal/mol)',80,20);ctx.fillText((hi-lo).toFixed(1),5,50);ctx.fillText('0',35,405)}document.getElementById('h').oninput=draw;document.getElementById('a').oninput=draw;draw();</script>'''
    template=template.replace('HVALUE',str(coupling)).replace('AVALUE',str(offset)).replace('DATA',json.dumps(rows))
    (output/'teaching_evb.html').write_text(template)
    print(json.dumps(summary,indent=2));return rows,summary
