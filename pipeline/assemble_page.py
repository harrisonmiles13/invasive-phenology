"""Assemble site/index.html and site/guides.html from their .src.html bodies.

The source is a page body (the form the Claude Artifact host expects). For
GitHub Pages it is wrapped in a full HTML document, with Leaflet's stylesheet
inlined and Leaflet's script loaded from cdnjs. `--inline` additionally
inlines the app scripts, for hosts that only serve a single page.

    python pipeline/assemble_page.py                 # -> site/index.html, site/guides.html
    python pipeline/assemble_page.py --inline OUT    # body-only, scripts inlined
"""
import os
import re
import sys

SITE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "site")
LEAFLET_JS = ('<script src="https://cdnjs.cloudflare.com/ajax/libs/leaflet/1.9.4/leaflet.js" '
              'crossorigin="anonymous" referrerpolicy="no-referrer"></script>')
APP_JS = ["js/model.js", "js/app.js"]


def read(rel):
    with open(os.path.join(SITE, rel), encoding="utf-8") as f:
        return f.read()


def body(inline, name="index"):
    if name != "index":
        return read(f"{name}.src.html")
    css = read("vendor/leaflet-1.9.4.css")
    css = re.sub(r"/\*.*?\*/", "", css, flags=re.S)
    css = re.sub(r"\n\s*\n", "\n", css)
    src = read("index.src.html").replace("/*LEAFLET_CSS*/", "/* Leaflet 1.9.4 (BSD-2-Clause) */\n" + css.strip())
    if inline:
        scripts = [LEAFLET_JS] + [f"<script>\n{read(p)}\n</script>" for p in APP_JS]
    else:
        scripts = [LEAFLET_JS] + [f'<script src="{p}"></script>' for p in APP_JS]
    return src.replace("<!--SCRIPTS-->", "\n".join(scripts))


def main():
    if len(sys.argv) == 3 and sys.argv[1] == "--inline":
        with open(sys.argv[2], "w", encoding="utf-8") as f:
            f.write(body(inline=True))
        return
    for name in ("index", "guides"):
        page = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
                '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n')
        head, rest = body(False, name).split("</style>", 1)
        page += head + "</style>\n</head>\n<body>\n" + rest.strip() + "\n</body>\n</html>\n"
        with open(os.path.join(SITE, f"{name}.html"), "w", encoding="utf-8") as f:
            f.write(page)
        print(f"wrote site/{name}.html")


if __name__ == "__main__":
    main()
