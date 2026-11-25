"""
render.py

Renders deliverable1.md to deliverable1.docx with pandoc, then sets the page
and type the brief asks for: US Letter, 1-inch margins, 11 pt sans body, single
spacing, compact headings, 9 pt captions, and a bordered 8 pt table.

Run from the repo root:

    python deliverable_1/render.py
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

D = Path(__file__).resolve().parent
SRC, OUT = D / "deliverable1.md", D / "deliverable1.docx"


def style_block(xml: str, sid: str) -> re.Match:
    m = re.search(rf'<w:style [^>]*w:styleId="{sid}".*?</w:style>', xml, flags=re.S)
    assert m, sid
    return m


def set_style(xml: str, sid: str, ppr: str | None = None, rpr: str | None = None) -> str:
    m = style_block(xml, sid)
    block = m.group()
    if ppr is not None:
        block = re.sub(r"<w:pPr>.*?</w:pPr>|<w:pPr\s*/>", "", block, flags=re.S)
        block = block.replace("</w:style>", f"<w:pPr>{ppr}</w:pPr></w:style>")
    if rpr is not None:
        block = re.sub(r"<w:rPr>.*?</w:rPr>|<w:rPr\s*/>", "", block, flags=re.S)
        block = block.replace("</w:style>", f"<w:rPr>{rpr}</w:rPr></w:style>")
    return xml[:m.start()] + block + xml[m.end():]


def fix_styles(s: str) -> str:
    font = '<w:rFonts w:ascii="Calibri" w:hAnsi="Calibri" w:eastAsia="Calibri" w:cs="Calibri"/>'
    s = re.sub(r"<w:rPrDefault>.*?</w:rPrDefault>",
               f'<w:rPrDefault><w:rPr>{font}<w:sz w:val="22"/><w:szCs w:val="22"/>'
               '<w:lang w:val="en-US"/></w:rPr></w:rPrDefault>', s, flags=re.S)
    s = re.sub(r"<w:pPrDefault>.*?</w:pPrDefault>",
               '<w:pPrDefault><w:pPr><w:spacing w:after="100" w:line="240" w:lineRule="auto"/>'
               '</w:pPr></w:pPrDefault>', s, flags=re.S)
    body = '<w:spacing w:before="0" w:after="100" w:line="240" w:lineRule="auto"/><w:jc w:val="both"/>'
    s = set_style(s, "BodyText", ppr=body)
    s = set_style(s, "Compact", ppr='<w:spacing w:before="20" w:after="20" w:line="240" w:lineRule="auto"/>')
    s = set_style(s, "Title", ppr='<w:spacing w:after="40"/><w:jc w:val="center"/>',
                  rpr=f'{font}<w:b/><w:color w:val="1F3864"/><w:sz w:val="30"/><w:szCs w:val="30"/>')
    for sid, sz in (("Subtitle", 22), ("Author", 22), ("Date", 20)):
        s = set_style(s, sid, ppr='<w:spacing w:after="20"/><w:jc w:val="center"/>',
                      rpr=f'{font}<w:color w:val="404040"/><w:sz w:val="{sz}"/><w:szCs w:val="{sz}"/>')
    s = set_style(s, "Heading1", ppr='<w:keepNext/><w:keepLines/><w:spacing w:before="200" w:after="60"/>'
                                     '<w:outlineLvl w:val="0"/>',
                  rpr=f'{font}<w:b/><w:color w:val="1F3864"/><w:sz w:val="26"/><w:szCs w:val="26"/>')
    s = set_style(s, "Caption", ppr='<w:spacing w:before="40" w:after="140" w:line="240" w:lineRule="auto"/>'
                                    '<w:jc w:val="both"/>',
                  rpr='<w:sz w:val="18"/><w:szCs w:val="18"/>')
    s = set_style(s, "TableCaption", ppr='<w:keepNext/><w:spacing w:before="120" w:after="60"/>',
                  rpr='<w:sz w:val="18"/><w:szCs w:val="18"/>')
    s = set_style(s, "CaptionedFigure", ppr='<w:keepNext/><w:spacing w:before="80" w:after="0"/>'
                                            '<w:jc w:val="center"/>')
    s = set_style(s, "Bibliography", ppr='<w:spacing w:after="60"/><w:ind w:left="360" w:hanging="360"/>',
                  rpr='<w:sz w:val="19"/><w:szCs w:val="19"/>')
    return s


def fix_document(x: str) -> str:
    borders = ('<w:tblBorders>' + "".join(f'<w:{e} w:val="single" w:sz="4" w:space="0" w:color="8C8C8C"/>'
                                          for e in ("top", "left", "bottom", "right", "insideH", "insideV"))
               + '</w:tblBorders>')

    def table(m):
        t = m.group()
        # schema order inside tblPr: ... tblW, jc, tblInd, tblBorders, shd, tblLayout, tblCellMar, tblLook
        pr = re.search(r"<w:tblPr>.*?</w:tblPr>", t, flags=re.S).group()
        nxt = re.search(r"<w:(shd|tblLayout|tblCellMar|tblLook)\b", pr)
        new_pr = pr[:nxt.start()] + borders + pr[nxt.start():] if nxt else pr.replace("</w:tblPr>", borders + "</w:tblPr>")
        t = t.replace(pr, new_pr, 1)
        t = re.sub(r"<w:tr>(?!<w:trPr>)", "<w:tr><w:trPr><w:cantSplit/></w:trPr>", t)
        t = t.replace("<w:trPr>", "<w:trPr><w:cantSplit/>").replace("<w:cantSplit/><w:cantSplit/>", "<w:cantSplit/>")
        t = re.sub(r'(<w:pStyle w:val="Compact"\s*/>)', r"\1<w:keepNext/>", t)
        t = t.replace("<w:rPr>", '<w:rPr><w:sz w:val="16"/><w:szCs w:val="16"/>')
        t = re.sub(r"<w:r>(?!<w:rPr>)", '<w:r><w:rPr><w:sz w:val="16"/><w:szCs w:val="16"/></w:rPr>', t)
        return t

    x = re.sub(r"<w:tbl>.*?</w:tbl>", table, x, flags=re.S)
    x = x.replace("</w:tbl>", '</w:tbl><w:p><w:pPr><w:spacing w:before="0" w:after="0" w:line="200" w:lineRule="exact"/></w:pPr></w:p>')
    # reference paragraphs (after the References heading) use a hanging indent
    head = x.find(">References<")
    if head > 0:
        tail = x[head:]
        tail = re.sub(r'<w:pStyle w:val="(BodyText|FirstParagraph)"\s*/>', '<w:pStyle w:val="Bibliography"/>', tail)
        tail = re.sub(r"<w:p>(?!<w:pPr>)", '<w:p><w:pPr><w:pStyle w:val="Bibliography"/></w:pPr>', tail)
        x = x[:head] + tail
    sect = ('<w:sectPr><w:pgSz w:w="12240" w:h="15840"/><w:pgMar w:top="1440" w:right="1440" '
            'w:bottom="1440" w:left="1440" w:header="708" w:footer="708" w:gutter="0"/></w:sectPr>')
    x = re.sub(r"<w:sectPr>.*?</w:sectPr>|<w:sectPr\s*/>", "", x, flags=re.S)
    return x.replace("</w:body>", sect + "</w:body>")


def main() -> None:
    subprocess.run(["pandoc", str(SRC), "-f", "markdown", "-t", "docx", "--resource-path", str(D),
                    "-o", str(OUT)], check=True)
    with tempfile.TemporaryDirectory() as tmp:
        new = Path(tmp) / "out.docx"
        with zipfile.ZipFile(OUT) as zi, zipfile.ZipFile(new, "w", zipfile.ZIP_DEFLATED) as zo:
            for item in zi.infolist():
                data = zi.read(item.filename)
                if item.filename == "word/styles.xml":
                    data = fix_styles(data.decode()).encode()
                elif item.filename == "word/document.xml":
                    data = fix_document(data.decode()).encode()
                elif item.filename == "[Content_Types].xml" and b'Extension="png"' not in data:
                    data = data.replace(b'<Default Extension="xml"',
                                        b'<Default Extension="png" ContentType="image/png" /><Default Extension="xml"')
                zo.writestr(item, data)
        shutil.move(new, OUT)
    print("wrote", OUT.name)


if __name__ == "__main__":
    main()
