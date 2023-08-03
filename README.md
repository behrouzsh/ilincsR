# ilincsR

### Building dev version
### New fast mario testing
docker build --no-cache --rm . -t shiny.ilincs.org:5000/ilincsr-dev-m:2022-11-14

### 2020
docker build --no-cache --rm . -t shiny.ilincs.org:5000/ilincsr-dev:2020-04-08

docker push shiny.ilincs.org:5000/ilincsr-dev:2020-04-08

####Running dev on gimm12

docker run -d -p 2009:80 -p 209:8006 -p 256:443 -e dport=3306 -e dhost=gimm2.ketl.uc.edu -e chost=dev.ilincs.org -e sigDB=ilincs_sigs  -v /mnt/raid/tmp/:/mnt/raid/tmp/ -v /opt/raid10/:/opt/raid10 --name ilincsr-dev shiny.ilincs.org:5000/ilincsr-dev:2020-04-08

#### Running dev on mario workstation

#### 2022 gimm14 New fast mario testing
docker run -d -p 8002:8006 -p 76:80 -p 439:443 -e dport=4040 -e dhost=eh3.ketl.uc.edu -e chost=eh332new.ketl.uc.edu -e sigDB=ilincs_new   -v /mnt/raid/tmp/:/mnt/raid/tmp/ -v /opt/raid10/:/opt/raid10 --name ilincsr-dev shiny.ilincs.org:5000/ilincsr-dev-m:2022-11-14

#### New fast mario testing
docker run -d -p 8002:8006 -p 76:80 -p 439:443 -e dport=4040 -e dhost=eh3.ketl.uc.edu -e chost=eh332new.ketl.uc.edu -e sigDB=ilincs_new   -v /mnt/raid/tmp/:/mnt/raid/tmp/ -v /opt/raid10/:/opt/raid10 --name ilincsr-dev shiny.ilincs.org:5000/ilincsr-dev-m:2022-11-14

#### 2020
docker run -d -p 8003:8006 -p 77:80 -p 440:443 -e dport=4040 -e dhost=eh3.ketl.uc.edu -e chost=eh332new.ketl.uc.edu -e sigDB=ilincs_new   -v /mnt/raid/tmp/:/mnt/raid/tmp/ -v /opt/raid10/:/opt/raid10 --name ilincsr-dev shiny.ilincs.org:5000/ilincsr-dev:2020-04-08

### To run for production

Go to Jenkins project BuildDocker_iLincsR2
- build container from the "master" branch
- click on the new build and promote to New
- click on the new build and promote to Prod


#Old - Running on mario workstation

docker run -d -p 8006:8006 -p 79:80 -p 443:443 -e dport=4040 -e dhost=eh3.ketl.uc.edu -e chost=eh332new.ketl.uc.edu -e sigDB=ilincs_new   -v /mnt/raid/tmp/:/mnt/raid/tmp/ -v /opt/raid10/:/opt/raid10 --name myilincsR shiny.ilincs.org:5000/ilincsr-opencpu:latest
