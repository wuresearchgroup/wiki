---
title: How to Use the Wiki
authors: Zhenghao Wu
comments: true
---

# How to Use the Wiki

The Wiki is written in markdown format. This wiki uses python-markdown as the markdown interpreter, supporting some markdown extension syntax. When editing markdown files locally, we recommend using [VSCode](https://code.visualstudio.com) or [Obsidian](https://obsidian.md/).

If you have any questions, you can provide feedback at <{{ config.repo_url }}/issues>.

> Sections marked with `*` in the document can be skipped.

## Have a question about a wiki page?

Please use the comment section at the bottom of the page and log in with your GitHub account to leave a comment. This feature is powered by [giscus](https://giscus.app/), which can automatically create a discussion thread for convenient interaction. Note: This function needs to be manually enabled by the contributor who created the page.

## Contribute to this website

Follow these steps to add or improve pages.

### 1) Fork and clone
- Fork the repository on GitHub.
- Clone your fork and enter the project directory:

```bash
git clone https://github.com/<your-username>/GroupWiki.git
cd GroupWiki
```

### 2) Create a working branch
```bash
git checkout -b docs/<topic-or-page>
```

### 3) Set up and preview locally
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
mkdocs serve
```
- Open the local preview URL shown in the terminal (usually `http://127.0.0.1:8000`).

### 4) Add or edit content
- Place new pages under `docs/wiki/` in an appropriate subfolder (e.g., `docs/wiki/new_comers/md/`).
- Start each page with front matter:

```yaml
---
title: Your Page Title
authors: Your Name
comments: true
---
```

- Images: put assets under `docs/images/` (keep a logical structure). Reference them with a path from the docs root, for example:
  - `![Caption](images/lammps/tutorial3/figures/PEG-distance.png)`

### 5) Add to the navigation
- Edit `mkdocs.yml` and add your page under `nav:` using paths relative to `docs/`, for example:

```yaml
nav:
  - Wiki:
    - How to Edit:
      - How to Use the Wiki: wiki/how_to_edit/howtodo.md
      - Your Page Title: wiki/<your-folder>/<your-file>.md
```

### 6) Commit and push
```bash
git add .
git commit -m "docs: add <your page title>"
git push origin docs/<topic-or-page>
```

### 7) Open a Pull Request
- Open a PR from your branch to the `main` branch of the upstream repository and briefly describe the changes.

### Content style checklist*
- Use clear, concise headings; prefer `##`/`###` levels.
- Use backticks for file names, commands, and code.
- Prefer relative links within the site and descriptive link text.
- Include alt text for images and keep images reasonably sized.
- Keep sections small and scannable; use lists and tables when helpful.