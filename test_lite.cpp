#include <cstdio>
#include <cstdlib>
#include <vector>
#include "lora_lite.h"
using namespace lora_lite;
static int g_count=0, g_crcok=0;
static void on_pkt(const Packet& p, void*){
    const char* st = p.status==PKT_DECODED?"OK ":p.status==PKT_CRC_FAIL?"BAD":"DET";
    if(p.status==PKT_DECODED) g_crcok++;
    if(p.status!=PKT_HEADER_ONLY) g_count++;
    printf("%s hdrlen=%2d len=%2d snr=%+.1f cfo=%+.1f sync=0x%02x  ", st, p.hdr_len, p.len, p.snr, p.cfo, p.sync);
    for(int i=0;i<p.len;i++) printf("%02x", p.data[i]);
    printf("\n");
}
int main(int argc, char** argv){
    if(argc<4){ fprintf(stderr,"usage: %s file.u8 sf os\n",argv[0]); return 1; }
    int sf=atoi(argv[2]);
    static Demod demod; demod.init(sf,1,true,0x12,on_pkt,nullptr);
    FILE* f=fopen(argv[1],"rb"); if(!f){perror("open");return 1;}
    std::vector<unsigned char> buf(1<<16); std::vector<Cx> samp; size_t nr;
    while((nr=fread(buf.data(),1,buf.size(),f))>0){
        int n=nr/2; samp.resize(n);
        for(int i=0;i<n;i++) samp[i]=Cx{((float)buf[2*i]-127.5f)/127.5f,((float)buf[2*i+1]-127.5f)/127.5f};
        demod.feed(samp.data(),n);
    }
    fclose(f);
    fprintf(stderr,"\nlora_lite: %d packets, %d CRC-valid\n", g_count, g_crcok);
    return 0;
}
