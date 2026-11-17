"""Read the complete native Q6 topology payload as a porting contract.

This is not yet a Q6-to-OpenMM Hamiltonian converter. Native units and index
spaces are preserved; solvent-boundary semantics require separate validation.
"""
from pathlib import Path
import math
from kemp_inputs import read,case_paths


class Reader:
    def __init__(self,text):self.lines=text.splitlines()
    def marker(self,text):
        matches=[i for i,line in enumerate(self.lines) if text in line]
        if len(matches)!=1:raise ValueError(f'Expected exactly one topology marker: {text}')
        i=matches[0];return i,self.lines[i]
    def counts(self,text):return [int(x) for x in self.marker(text)[1].split('=')[0].split()]
    def values(self,text,n,cast=float):
        i,_=self.marker(text);result=[]
        while len(result)<n:
            i+=1
            if i>=len(self.lines) or '=' in self.lines[i]:raise ValueError(f'Truncated topology block: {text}')
            result.extend(cast(x) for x in self.lines[i].split())
        if len(result)!=n:raise ValueError(f'Unexpected payload length: {text}')
        return result
    def rows(self,text,n,width):
        values=self.values(text,n*width,int)
        return [values[i:i+width] for i in range(0,len(values),width)]
    def parameter_rows(self,text):
        count=self.counts(text)[0];i,_=self.marker(text);result={}
        for line in self.lines[i+1:i+1+count]:
            fields=line.split();code=int(fields[0])
            if code in result:raise ValueError('Duplicate topology parameter code.')
            result[code]=list(map(float,fields[1:]))
        if len(result)!=count:raise ValueError('Incomplete parameter table.')
        return result
    def text_values(self,text,n,allowed):
        i,_=self.marker(text);characters=[]
        while len(characters)<n:
            i+=1
            if i>=len(self.lines) or '=' in self.lines[i]:raise ValueError(f'Truncated bitfield: {text}')
            characters.extend(self.lines[i].strip())
        if len(characters)!=n or any(c not in allowed for c in characters):raise ValueError('Malformed topology bitfield.')
        return characters


def parse_topology(text):
    r=Reader(text);n,solute,water_size,dielectric=r.counts('Total no. of atoms')
    # Numeric parameter-table and textual type-name markers differ only by colon.
    masses_line,_=r.marker('Masses:')
    type_header=r.lines[masses_line-3]
    ntypes=int(type_header.split('=')[0])
    coordinates=r.values('Total no. of atoms',3*n)
    iac=r.values("No. of integer atom codes",n,int)
    charges=r.values('No. of atomic charges',n)
    masses=r.values('Masses:',ntypes)
    nonbonded={label:r.values(label,ntypes) for label in ('sqrt (Aii) normal:','sqrt (Bii) normal:','sqrt (Aii) 1-4:','sqrt (Bii) 1-4:')}
    blocks={}
    for name,marker,width,pm in [('bonds','No. of bonds, no. of solute bonds',3,'No. of bond codes.'),('angles','No. of angles, no. of solute angles',4,'No. of angle codes.'),('torsions','No. of torsions, solute torsions',5,'No. of torsion codes.'),('impropers','No. of impropers, solute impr.',5,'No. of improper codes, type')]:
        count=r.counts(marker)[0];rows=r.rows(marker,count,width);parameters=r.parameter_rows(pm)
        for row in rows:
            if any(i<1 or i>n for i in row[:-1]) or row[-1] not in parameters:raise ValueError(f'Invalid {name} atom/type reference.')
        blocks[name]={'terms':rows,'parameter_types':parameters}
    # Q6 topo.f90 f04f99e: column-major temp_list(25,nat_solute),
    # list(delta, lower_index) represents pair(lower,lower+delta).
    pairs={}
    for name,marker,longmarker in [('one_four','No. of 1-4 neighbours.','No. of long 1-4 nbrs'),('excluded','No. of exclusions.','No. of long exclusions')]:
        declared=r.counts(marker)[0];bits=r.text_values(marker,25*solute,{'0','1'});short=[]
        for lower in range(1,solute+1):
            for delta in range(1,26):
                if bits[(lower-1)*25+delta-1]=='1':
                    if lower+delta>n:raise ValueError('Out-of-range native pair bit.')
                    short.append([lower,lower+delta])
        long_count=r.counts(longmarker)[0];long=r.rows(longmarker,long_count,2)
        if name=='one_four' and len(short)!=declared:raise ValueError('Decoded 1-4 count differs from native header.')
        # Qprep counts every bonded/angle event, including solvent/duplicates,
        # but serializes packed bits only for the solute. Preserve that distinction.
        if name=='excluded':
            for oxygen in range(solute+1,n+1,water_size):
                if water_size!=3:raise ValueError('Only supplied three-site waters are supported.')
                short.extend([[oxygen,oxygen+1],[oxygen,oxygen+2],[oxygen+1,oxygen+2]])
        pairs[name]={'short_pairs':short,'long_pairs':long,'native_header_event_count':declared}
    if any(i<1 or i>ntypes for i in iac):raise ValueError('Invalid atom type code.')
    if not all(math.isfinite(x) for x in (*coordinates,*charges,*masses)):raise ValueError('Nonfinite native input.')
    scale_line=r.marker('Electrostatic 1-4 scaling factor')[1];scaling,coulomb=map(float,scale_line.split('=')[0].split())
    excluded_bits=r.text_values('No. of excluded atoms (incl. water)',n,{'T','F'})
    center_line=r.marker('Solvent centre')[1];center=list(map(float,center_line.split('=')[0].split()))
    radii=list(map(float,r.marker('Exclusion, solvent radii')[1].split('=')[0].split()))
    return {'excluded_atom_indices_1based':[i+1 for i,x in enumerate(excluded_bits) if x=='T'],
            'solvent_center_angstrom':center,'solvent_radius_angstrom':radii[1],
            'schema':'q6_native_topology_port_contract_v1','atoms':n,'solute_atoms':solute,'water_size':water_size,
            'coordinates_angstrom':[coordinates[i:i+3] for i in range(0,len(coordinates),3)],'atom_type_codes':iac,
            'charges_e':charges,'type_masses_da':masses,'native_LJ_tables':nonbonded,'bonded':blocks,'pairs':pairs,
            'vdw_combination_rule':r.counts('vdW combination rule')[0], 'coulomb_14_scale':scaling,'coulomb_kcal_angstrom_mol_e2':coulomb,
            'solvent_dielectric':dielectric/1000,'converted_to_openmm':False,
            'note':'Complete native coefficient/index payload; no force generation, FEP patching, boundary implementation or scientific approval.'}


def audit_case_topology(case):
    data=parse_topology(read(case_paths(case)[2]))
    return {'case':case,'atoms':data['atoms'],'solute_atoms':data['solute_atoms'],
            'bonded_counts':{k:len(v['terms']) for k,v in data['bonded'].items()},
            'pair_counts':{k:len(v['short_pairs'])+len(v['long_pairs']) for k,v in data['pairs'].items()},
            'net_topology_charge_e':sum(data['charges_e']),'coulomb_native':data['coulomb_kcal_angstrom_mol_e2'],
            'vdw_combination_rule':data['vdw_combination_rule'],'converted_to_openmm':False}
