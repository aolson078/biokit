from flask import Flask, request, render_template
from werkzeug.utils import secure_filename
import os
from generate_tree import generate_tree

app = Flask(__name__)


@app.route('/')
def home():
    return render_template('index.html')


@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return 'File not uploaded'
    file = request.files['file']

    if file.filename == '':
        return 'No selected file'
    if file:
        filename = secure_filename(file.filename)
        filepath = os.path.join('uploads', filename)
        file.save(filepath)
        generate_tree(filepath)  # Ensure generate_tree uses the uploaded file
        return 'File uploaded and tree generated'

    return '0'


@app.route('/tree')
def display_tree():
    # Write logic for displaying tree

    return '0'


if __name__ == '__main__':
    app.run()
