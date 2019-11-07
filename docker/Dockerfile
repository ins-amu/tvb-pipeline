FROM gitlab.thevirtualbrain.org:5000/docker/pipeline/master
MAINTAINER Marmaduke Woodman <marmaduke.woodman@univ-amu.fr>

ADD ./ /opt/pipeline
WORKDIR /opt/pipeline
ENTRYPOINT ["make"]