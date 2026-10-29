"""Read-only retrieval of original papers for figure-layout inspection."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json
import requests

ROOT = Path(__file__).resolve().parent
PAPERS = {
    'alphafold3': 'https://www.nature.com/articles/s41586-024-07487-w.pdf',
    'rmsfnet': 'https://www.nature.com/articles/s41467-024-49858-x.pdf',
    'bioemu_accepted': 'https://www.microsoft.com/en-us/research/wp-content/uploads/2024/12/accepted_paper.pdf',
    'alphamissense_vor_mirror': 'https://static1.squarespace.com/static/603ab50ab81d5532a0a4a42b/t/651480e7c4cab752c53933f6/1695842536657/Cheng%2Bet%2Bal.pdf',
}

def retrieve(item):
    name, url = item
    try:
        response = requests.get(url, timeout=25)
        result = dict(name=name, url=url, status=response.status_code,
                      bytes=len(response.content), final_url=response.url)
        if response.status_code == 200 and response.content.startswith(b'%PDF'):
            target = ROOT / (name + '.pdf')
            target.write_bytes(response.content)
            result['local_path'] = str(target)
        return result
    except Exception as exc:
        return dict(name=name, url=url, error=str(exc))

if __name__ == '__main__':
    results = list(ThreadPoolExecutor(max_workers=3).map(retrieve, PAPERS.items()))
    (ROOT / 'paper_retrieval.json').write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(results, ensure_ascii=False, indent=2))
