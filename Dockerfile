# FracHinge web application (Streamlit). Works on Hugging Face Spaces (Docker SDK),
# Render, Railway, Google Cloud Run or any server with Docker.
FROM python:3.11-slim
WORKDIR /app
COPY requirements-app.txt .
RUN pip install --no-cache-dir -r requirements-app.txt
COPY frachinge ./frachinge
COPY app ./app
COPY examples ./examples
COPY docs/tables ./docs/tables
COPY .streamlit ./.streamlit
ENV PORT=8501
EXPOSE 8501
CMD streamlit run app/streamlit_app.py --server.port=${PORT} --server.address=0.0.0.0
