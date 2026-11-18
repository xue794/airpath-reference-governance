from pathlib import Path
import os,sys,shutil,json
H=Path(__file__).resolve().parent
skill=Path(r'C:\Users\IMUT\.codex\plugins\cache\openai-primary-runtime\documents\26.1104.11149\skills\documents')
os.environ['PATH']=r'C:\Users\IMUT\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin'+os.pathsep+os.environ['PATH']
sys.path.insert(0,str(skill));import render_docx as rd
def conversion(doc_path,user_profile,convert_tmp_dir,stem,verbose=False):
 target=Path(convert_tmp_dir)/(stem+'.pdf');shutil.copy2(H/'qa/word_export.pdf',target)
 return str(target),'Native Word read-only export; bundled rasterizer; bundled LibreOffice unavailable on host'
rd.convert_to_pdf=conversion
paths=rd.rasterize(str(H/'Aerospace_WP3反馈裕度传递核查_v1_0.docx'),str(H/'qa/report'),135,False,False)
(H/'qa/render_summary.json').write_text(json.dumps(dict(backend='Word read-only export + bundled render_docx.py rasterize',pages=len(paths),paths=paths),ensure_ascii=False,indent=2),encoding='utf-8')
print('RENDERED_PAGES',len(paths))
