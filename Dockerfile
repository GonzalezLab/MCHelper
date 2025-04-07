FROM dfam/tetools:latest AS builder

RUN apt update && apt upgrade -y

# Add software to path



# MCHelper installation

RUN apt install -y git unzip python3-pandas python3-opencv python3-psutil python3-sklearn r-base python3-opencv hmmer emboss
RUN pip install --break-system-packages pdf2image Bio cialign
