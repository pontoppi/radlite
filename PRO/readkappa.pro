function readkappa,file
openr,1,file
iformat=0
readf,1,iformat
case iformat of
    1: begin
        nf=0
        readf,1,nf
        data = dblarr(2,nf)
        readf,1,data
        lambda    = transpose(data[0,*])
        kappa_abs = transpose(data[1,*])
        kappa_sca = dblarr(nf)
        g_sca     = dblarr(nf)
    end
    2: begin
        nf=0
        readf,1,nf
        data = dblarr(3,nf)
        readf,1,data
        lambda    = transpose(data[0,*])
        kappa_abs = transpose(data[1,*])
        kappa_sca = transpose(data[2,*])
        g_sca     = dblarr(nf)
    end
    3: begin
        nf=0
        readf,1,nf
        data = dblarr(4,nf)
        readf,1,data
        lambda    = transpose(data[0,*])
        kappa_abs = transpose(data[1,*])
        kappa_sca = transpose(data[2,*])
        g_sca     = transpose(data[3,*])
    end
endcase
close,1	
freq = 2.9979d14/lambda
return,{freq:freq,lambda:lambda,kappa_abs:kappa_abs,kappa_sca:kappa_sca,$
           g_sca:g_sca}
end
