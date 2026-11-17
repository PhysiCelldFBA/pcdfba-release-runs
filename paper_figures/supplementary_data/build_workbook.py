#!/usr/bin/env python3
"""Build Supplementary_Data.xlsx from the tables of 06_supporting_material.tex.

Usage:  python build_workbook.py [TEX [OUT]]
  TEX: LaTeX source of the tables (default tables_source.tex, extracted from the supplement)
  OUT: output workbook (default Supplementary_Data.xlsx)
Requires openpyxl.
Tables with where == 'workbook' or 'both' in table_map.py become sheets.  A README sheet lists them.
"""
import re, sys
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter
from table_map import TABLES

SRC = sys.argv[1] if len(sys.argv) > 1 else 'tables_source.tex'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'Supplementary_Data.xlsx'

def group(s, i):
    """text inside the brace group starting at s[i]=='{' -> (inner, index after)"""
    d = 0
    for j in range(i, len(s)):
        if s[j] == '{' and s[j-1] != '\\': d += 1
        elif s[j] == '}' and s[j-1] != '\\':
            d -= 1
            if d == 0: return s[i+1:j], j+1
    raise ValueError('unbalanced')

def unwrap(s, cmd, nopt=0, keep_arg=-1):
    """replace \\cmd[opt]{a}{b}... by the argument number keep_arg (0-based); keep_arg=-1: last"""
    pat = re.compile(r'\\' + cmd + r'(?![A-Za-z])')
    while True:
        m = pat.search(s)
        if not m: return s
        i = m.end()
        if i < len(s) and s[i] == '[':
            i = s.index(']', i) + 1
        args = []
        for _ in range(nopt):
            while s[i] == ' ': i += 1
            a, i = group(s, i); args.append(a)
        s = s[:m.start()] + (args[keep_arg] if args else '') + s[i:]

SYM = {r'\Delta':'Δ', r'\mu':'µ', r'\times':'×', r'\le':'≤', r'\leq':'≤', r'\ge':'≥', r'\geq':'≥',
       r'\approx':'≈', r'\pm':'±', r'\sim':'~', r'\to':'→', r'\rightarrow':'→', r'\alpha':'α', r'\beta':'β',
       r'\max':'max', r'\min':'min', r'\%':'%', r'\_':'_', r'\&':'&', r'\#':'#', r'\,':'', r'\;':' ', r'\ ':' ',
       r'\cdot':'·', r'\mathrm':'', r'\textit':'', r'\textbf':'', r'\texttt':'', r'\text':'', r'\mathit':''}

def clean(s):
    s = re.sub(r'\\chcolor\b', '', s)
    s = re.sub(r'\\parbox\[[^\]]*\]\{[^}]*\}', '', s)
    s = re.sub(r'\\(Delta|mu|alpha|beta)(?![A-Za-z])\s*', lambda m: {'Delta':'Δ','mu':'µ','alpha':'α','beta':'β'}[m.group(1)], s)
    for c, n in [('chadded', 1), ('chdeleted', 1), ('chreplaced', 2), ('chhighlight', 1)]:
        pass
    # tracked-change wrappers: keep the new text
    s = unwrap(s, 'chadded', 1)
    s = unwrap(s, 'chreplaced', 2, 0)
    s = unwrap(s, 'chdeleted', 1) if False else re.sub(r'\\chdeleted(\[[^\]]*\])?\{[^{}]*\}', '', s)
    s = re.sub(r'\\(textbf|textit|mathrm|texttt|text|mathit|emph)\s*\{', '{', s)
    for k in sorted(SYM, key=len, reverse=True):
        s = s.replace(k, SYM[k])
    s = re.sub(r'\$\^\{?([^$}]*)\}?\$', r'^\1', s)       # simple superscripts
    s = s.replace('$', '')
    s = re.sub(r'_\{([^}]*)\}', r'_\1', s)
    s = re.sub(r'\^\{([^}]*)\}', r'^\1', s)
    s = re.sub(r'(?<=[A-Za-z])_(\d)', lambda m: '₀₁₂₃₄₅₆₇₈₉'[int(m.group(1))], s)
    s = s.replace('~', ' ').replace('---', '—').replace('--', '–')
    s = re.sub(r'\\[a-zA-Z]+\*?', '', s)
    s = s.replace('{', '').replace('}', '').replace('\\', '')
    return re.sub(r'\s+', ' ', s).strip()

def split_top(s, sep):
    out, d, cur, i = [], 0, '', 0
    while i < len(s):
        c = s[i]
        if c == '\\' and s[i:i+2] == '\\\\' and sep == '\\\\' and d == 0:
            out.append(cur); cur = ''; i += 2; continue
        if c == '\\' and i+1 < len(s):
            cur += s[i:i+2]; i += 2; continue
        if c == '{': d += 1
        if c == '}': d -= 1
        if c == '&' and sep == '&' and d == 0:
            out.append(cur); cur = ''; i += 1; continue
        cur += c; i += 1
    out.append(cur)
    return out

def parse_table(env):
    cap = ''
    m = re.search(r'\\caption\s*\{', env)
    if m: cap, _ = group(env, m.end()-1)
    b = env.index('\\begin{tabular}'); e = env.index('\\end{tabular}')
    body = env[b:e]
    body = body[body.index('}', body.index('\\begin{tabular}{') + len('\\begin{tabular}')) + 1:] if False else body
    # skip the column spec (balanced braces)
    k = body.index('{', len('\\begin{tabular}')); _, k2 = group(body, k); body = body[k2:]
    body = re.sub(r'\\(toprule|bottomrule)', '', body)
    body = re.sub(r'\\cmidrule(\([^)]*\))?\{[^}]*\}', '', body)
    parts = re.split(r'(\\midrule)', body)
    rows = []; header_done = False; nhead = 0
    for p in parts:
        if p == '\\midrule':
            header_done = True; continue
        for r in split_top(p, '\\\\'):
            if not r.strip(): continue
            cells = []
            for c in split_top(r, '&'):
                c = c.strip(); span = 1
                mm = re.match(r'\\multicolumn\{(\d+)\}\{[^}]*\}\s*\{', c)
                if mm:
                    span = int(mm.group(1)); c, _ = group(c, mm.end()-1)
                bold = '\\textbf' in c
                cells.append((clean(c), span, bold))
            rows.append((cells, not header_done))
    return clean(cap), rows

def to_num(t):
    if re.fullmatch(r'-?\d{1,3}(,\d{3})+', t): return int(t.replace(',', ''))
    if re.fullmatch(r'-?\d+', t): return int(t)
    if re.fullmatch(r'-?\d*\.\d+', t): return float(t)
    return t

SECTION = {'UnitTest_design':'S4','UnitTest_volume_deviation':'S4','UnitTest_voxel_vs_cell':'S4','Tumour_dt_dx':'S7',
           'Consortium_dt_dx':'S7','Consortium_nocaps_dt':'S7','Motility_dt_dx':'S7','Cap_number':'S7',
           'Tumour_supply_grid':'S10','Consortium_sensitivity':'S10','Cost_scaling':'S11','Cost_24h_estimate':'S11'}

def main():
    src = open(SRC).read()
    wb = Workbook(); idx = wb.active; idx.title = 'Index'
    index = []
    for old, label, where, sheet in TABLES:
        if where == 'pdf': continue
        i = src.index('\\label{tab:%s}' % label)
        a = src.rfind('\\begin{table', 0, i); b = src.index('\\end{table', i)
        env = src[a:b]
        cap, rows = parse_table(env)
        ws = wb.create_sheet(sheet)
        ws.append([cap]); ws['A1'].font = Font(bold=True)
        ws.append([])
        r0 = 3
        ncol = 1
        for cells, hdr in rows:
            col = 1
            for t, span, bold in cells:
                c = ws.cell(row=r0, column=col, value=to_num(t))
                c.alignment = Alignment(wrap_text=True, vertical='top')
                if hdr or bold: c.font = Font(bold=True)
                if span > 1: ws.merge_cells(start_row=r0, start_column=col, end_row=r0, end_column=col+span-1)
                col += span
            ncol = max(ncol, col-1); r0 += 1
        for j in range(1, ncol+1):
            w = max((len(str(ws.cell(row=r, column=j).value or '')) for r in range(3, r0)), default=10)
            ws.column_dimensions[get_column_letter(j)].width = min(max(12, w+2), 48)
        ws.freeze_panes = 'A4'
        r = re.search(r'\\ref\{tab:%s\}' % label, src)
        sec = ''
        if r:
            hs = re.findall(r'\\subsection\*\{(?:\\chadded\[id=mpdl\]\{)?(S\d+(?:\.\d+)?)\.', src[:r.start()])
            sec = hs[-1] if hs else ''
        index.append((sheet, cap, SECTION[sheet]))
    idx.append(['Sheet', 'Contents', 'Supplementary Methods section']); 
    for c in idx[1]: c.font = Font(bold=True)
    for s_, c, sec in index: idx.append([s_, c, sec])
    idx.column_dimensions['A'].width = 30; idx.column_dimensions['B'].width = 120; idx.column_dimensions['C'].width = 18
    for row in idx.iter_rows(min_row=2):
        row[1].alignment = Alignment(wrap_text=True, vertical='top')
    wb.save(OUT)
    print('wrote', OUT, [x[0] for x in index])

if __name__ == '__main__':
    main()
