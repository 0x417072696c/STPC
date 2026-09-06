from flask import Flask, request, jsonify, render_template_string
import allymind_core, cascade, pdb_parser
import os

app = Flask(__name__)

HTML_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <title>AllyMind – Mutation Analysis</title>
    <meta charset="UTF-8">
    <style>
        body { font-family: 'Segoe UI', sans-serif; background: #1a1a1a; color: #eff0f1; padding: 40px; }
        h1 { color: #3daee9; text-align: center; }
        form { max-width: 500px; margin: 30px auto; background: #2d2d2d; padding: 25px; border-radius: 10px; }
        label { display: block; margin-top: 15px; color: #888; }
        input, button { width: 100%; padding: 10px; margin-top: 5px; border-radius: 5px; border: 1px solid #3daee9; background: #1a1a1a; color: #eff0f1; }
        button { background: #3daee9; color: white; font-weight: bold; cursor: pointer; margin-top: 20px; }
        button:hover { background: #2980b9; }
        #result { max-width: 600px; margin: 30px auto; background: #2d2d2d; padding: 20px; border-radius: 10px; white-space: pre-wrap; display: none; }
        .footer { text-align: center; margin-top: 40px; font-size: 0.8rem; color: #555; }
    </style>
</head>
<body>
    <h1>AllyMind Web Server</h1>
    <form id="mutationForm">
        <label>UniProt ID or gene symbol:</label>
        <input type="text" id="uniprot" value="P00441" required>
        <label>Mutation (e.g., A4V):</label>
        <input type="text" id="mutation" value="A4V" required>
        <label>PDB ID (optional):</label>
        <input type="text" id="pdb" value="2C9V">
        <button type="submit">Analyze</button>
    </form>
    <div id="result"></div>
    <div class="footer">
        AllyMind Web Server – Confidential. Valentin Lebedkin, 2026.
    </div>
    <script>
        document.getElementById('mutationForm').addEventListener('submit', async (e) => {
            e.preventDefault();
            const uniprot = document.getElementById('uniprot').value;
            const mutation = document.getElementById('mutation').value;
            const pdb = document.getElementById('pdb').value;
            const resultDiv = document.getElementById('result');
            resultDiv.style.display = 'block';
            resultDiv.textContent = 'Calculating...';
            try {
                const response = await fetch('/analyze', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ uniprot, mutation, pdb })
                });
                const data = await response.json();
                if (data.error) {
                    resultDiv.textContent = 'Error: ' + data.error;
                } else {
                    resultDiv.textContent = JSON.stringify(data, null, 2);
                }
            } catch (err) {
                resultDiv.textContent = 'Connection error: ' + err.message;
            }
        });
    </script>
</body>
</html>
"""

@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)

@app.route('/analyze', methods=['POST'])
def analyze():
    data = request.get_json()
    uniprot = data.get('uniprot', '').strip().upper()
    mutation = data.get('mutation', '').strip().upper()
    pdb_id = data.get('pdb', '').strip()

    if not uniprot or not mutation:
        return jsonify({'error': 'UniProt ID and mutation are required'}), 400

    seq, uid = allymind_core.load_sequence(uniprot)
    if not seq:
        return jsonify({'error': f'Protein {uniprot} not found'}), 404

    ss_map = None
    if pdb_id:
        pdb_file = pdb_parser.download_pdb(pdb_id)
        if pdb_file:
            ss_map = pdb_parser.extract_aligned_ss_map(seq, pdb_file, chain='A')

    try:
        report = allymind_core.analyze_mutation(seq, mutation, ss_map=ss_map)
    except Exception as e:
        return jsonify({'error': f'Analysis failed: {str(e)}'}), 500

    wt = mutation[0]
    pos = int(mutation[1:-1]) - 1
    alt = mutation[-1]
    depth = report.get('depth', 0)
    ddG = report.get('kinetics', {}).get('delta_G', 0) or 0
    thr = report.get('therapeutic_threshold', 0) or 0
    ss = report.get('secondary_structure', 'loop') or 'loop'
    exposed = cascade.is_trap_exposed(seq, pos, wt, alt, ss_map=ss_map)
    classification = "open" if exposed else "hidden"

    if classification == "open":
        strategy = "Small molecule" if depth < 1.0 else "Chaperone"
    else:
        strategy = "Monitor" if depth < thr else "Chaperone (hidden)"

    return jsonify({
        'protein': uid,
        'mutation': mutation,
        'depth': round(depth, 4),
        'ddG': round(ddG, 4),
        'threshold_J': round(thr, 4),
        'structure': ss,
        'classification': classification,
        'strategy': strategy
    })

if __name__ == '__main__':
    print("AllyMind Web Server starting...")
    print("Open http://127.0.0.1:5000 in your browser")
    app.run(host='0.0.0.0', port=5000, debug=True)
