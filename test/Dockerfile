FROM python:latest
WORKDIR /usr/app/src
COPY . .
RUN pip3 install --upgrade pip --no-cache-dir \
    && pip3 install -r ./requirements.txt --no-cache-dir 
CMD [ "python", "./delete_vm.py"]
