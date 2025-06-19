# Entry point for running the upgraded application
from flask_bio_app import create_app

def deploy():
    app = create_app()
    app.run()

if __name__ == '__main__':
    deploy()
