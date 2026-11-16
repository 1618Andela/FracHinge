# Publishing the FracHinge web application / Publicar la aplicación web

The app (`app/streamlit_app.py`) runs in any browser once deployed. Three free or low-cost routes:

## A. Streamlit Community Cloud (simplest, free)
1. Create a public GitHub repository and push this folder (`git init; git add .; git commit -m "FracHinge 1.0.1"; git push`).
2. Go to https://share.streamlit.io, sign in with GitHub, click **Create app**.
3. Repository: `<user>/FracHinge`, branch `main`, main file path: `app/streamlit_app.py`.
4. Deploy. The public URL is `https://<name>.streamlit.app`. Dependencies are read from `requirements.txt`.

## B. Hugging Face Spaces (free, Docker)
1. Create a Space with SDK **Docker**.
2. Push this repository to the Space (the `Dockerfile` is ready; port 8501 must be declared
   by adding `app_port: 8501` to the Space README front matter).

## C. Own server / institutional server
```
docker build -t frachinge .
docker run -p 8501:8501 frachinge
```
Then open http://<server>:8501 (put it behind the university's HTTPS reverse proxy).

## Before going public
- Add the Zenodo DOI and the paper reference to the sidebar and README.
- Keep the disclaimer: the tool supports, but does not replace, the designer's verification against EN 1992-1-1:2023.

---

# Publicar la aplicación (resumen en español)
- **Streamlit Community Cloud**: subir el repositorio a GitHub, entrar en share.streamlit.io, *Create app*, fichero principal `app/streamlit_app.py`. URL pública `https://<nombre>.streamlit.app`.
- **Hugging Face Spaces**: Space de tipo Docker con el `Dockerfile` incluido (`app_port: 8501`).
- **Servidor propio**: `docker build -t frachinge . && docker run -p 8501:8501 frachinge`.
