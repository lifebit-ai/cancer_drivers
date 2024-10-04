################## BASE IMAGE ######################
FROM ensemblorg/ensembl-vep@sha256:70d64039e7b575801390296fa77bd78b983c94ab204baaa656425d8d4b738dac

################## METADATA ######################
LABEL software="vep" \
version="110.1" \
software.version="v110.1" \
about.home="https://www.ensembl.org/info/docs/tools/vep/" \
maintainer="Eugene Nadezhdin <eugene@lifebit.ai>"

################## INSTALLATION ######################
USER root

#Installing miniconda to install Samtools
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
apt-get install -y \
wget \
bzip2 \
ca-certificates \
curl \
git \
procps \
#Install dependencies for vep plugins
pkg-config \
libgd-dev \
lib32stdc++6 \
sqlite3 && \
# pkg-config libgd-dev are needed for GD perl module.
# lib32stdc++6 is required to run the binary of GeneSplicer plugin, that was built under older Linux version.
# sqlite3 is required for vep plugin LocalID
apt-get clean && \
rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh -O ~/miniconda.sh && \
/bin/bash ~/miniconda.sh -b -p /opt/conda && \
rm ~/miniconda.sh && \
/opt/conda/bin/conda clean -tipsy && \
ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
echo "conda activate base" >> ~/.bashrc


COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/cromwell-env/bin:$PATH

#Install dependencies for vep plugins
#RUN apt-get update && apt-get install -y pkg-config libgd-dev lib32stdc++6 sqlite3
# pkg-config libgd-dev are needed for GD perl module.
# lib32stdc++6 is required to run the binary of GeneSplicer plugin, that was built under older Linux version.
# sqlite3 is required for vep plugin LocalID
RUN cpanm install Math::CDF Net::FTP List::MoreUtils POSIX GD::Polygon GD DBI Text::CSV Scalar::Util FileHandle Cwd Digest::MD5 File::Basename LWP::UserAgent List::Util Data::Dumper Storable DBD::SQLite
# DBD::SQLite is needed for vep plugin LocalID

# vep module PON_P2 requires python module suds, which in term requires module client
#RUN pip install client
RUN pip install suds-py3

USER root

WORKDIR /data/

CMD ["bash"]

Close
