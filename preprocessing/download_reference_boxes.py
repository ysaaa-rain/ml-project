"""Download official RegulonDB short-box/citation snapshot locally only."""
from pathlib import Path
import subprocess,json,hashlib
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'tmp/reference_library_v1_20261008';URL='https://regulondb.ccg.unam.mx/graphql'
def query(q):
 result=json.loads(subprocess.check_output(['curl','--fail','-sS','--max-time','60',URL,'-H','Content-Type: application/json','--data-binary',json.dumps({'query':q})]));assert 'errors' not in result,result;return result
CITATIONS='citations{evidence{code name type} publication{pmid title citation url year}}'
QUERY='''{getAllOperon(limit:10000,page:0){pagination{currentPage lastPage totalResults hasNextPage} data{_id operon{_id name strand} transcriptionUnits{_id name promoter{_id name confidenceLevel note sequence boxes{leftEndPosition rightEndPosition sequence type} transcriptionStartSite{leftEndPosition rightEndPosition} CITATIONS bindsSigmaFactor{_id name abbreviatedName CITATIONS}}}}}}'''.replace('CITATIONS',CITATIONS)
def main():
 BASE.mkdir(parents=True,exist_ok=True)
 for p in ['regulondb_operons_all.json','regulondb_promoters.json','regulondb_version.json']:assert not (BASE/p).exists(),'Frozen snapshot exists; do not overwrite'
 version=query('{getDatabaseInfo{regulonDBVersion ecocycVersion genomeVersion releaseDate}}');data=query(QUERY);d=data['data']['getAllOperon'];assert not d['pagination']['hasNextPage'] and len(d['data'])==d['pagination']['totalResults']
 records={}
 for o in d['data']:
  for tu in o.get('transcriptionUnits') or []:
   if tu.get('promoter'):
    r=tu['promoter'];r['operon_strand']=o['operon']['strand'];records[r['_id']]=r
 for name,value in [('regulondb_operons_all.json',data),('regulondb_promoters.json',list(records.values())),('regulondb_version.json',version)]:
  (BASE/name).write_text(json.dumps(value,indent=2))
 m={'endpoint':URL,'query':QUERY,'access_date':'2026-10-08','operons':len(d['data']),'operon_associated_promoters':len(records),'does_not_certify_orphan_promoter_coverage':True,'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in BASE.glob('regulondb*.json')}};(BASE/'download_manifest.json').write_text(json.dumps(m,indent=2)+'\n')
if __name__=='__main__':main()
