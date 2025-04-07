FROM dfam/tetools:latest AS builder

RUN apt update && apt upgrade -y

# Add software to path



# MCHelper installation
RUN apt install -y git unzip python3-pandas python3-opencv python3-psutil python3-sklearn r-base python3-opencv hmmer emboss
RUN pip install --break-system-packages pdf2image Bio cialign
# Add dependencies to path
ENV PATH="$PATH:/opt/cd-hit:/opt/mafft/bin"

RUN git clone https://github.com/GonzalezLab/MCHelper /opt/mchelper
WORKDIR /opt/mchelper
RUN cd db && unzip '*.zip'
RUN bash
