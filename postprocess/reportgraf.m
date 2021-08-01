function reportgraf( caso, malha )
%
for ds = [1 3]
    
    if ds==1, metsat='SEQ'; elseif ds==2, metsat='Referência'; elseif ds==3, metsat='Streamline'; end 
    
    formethod = 'MPFAD';
    
    if exist(sprintf('%s\\%s\\%s\\%s\\Results\\Cumulateoil.mat',caso,malha,formethod,metsat),'file') ~= 0
        load(sprintf('%s\\%s\\%s\\%s\\Results\\Cumulateoil',caso,malha,formethod,metsat));
    end

    if exist(sprintf('%s\\%s\\%s\\%s\\Results\\Watercut.mat',caso,malha,formethod,metsat),'file') ~= 0
        load(sprintf('%s\\%s\\%s\\%s\\Results\\Watercut',caso,malha,formethod,metsat));
    end
    
    if exist(sprintf('%s\\%s\\%s\\%s\\Results\\VPI.mat',caso,malha,formethod,metsat),'file') ~= 0
        load(sprintf('%s\\%s\\%s\\%s\\Results\\VPI',caso,malha,formethod,metsat));
    end
    
    if ds==1
        cumseq = cumulateoil;
        wcutseq = watercut;
        vpiseq = VPI;
    elseif ds==3
        cumsl = cumulateoil;
        wcutsl = watercut;
        vpisl = VPI;
    elseif ds==2
        cumbr = CumOil;
        wcutbr = Watercut;
        vpibr = VPI;
    end
    
end

figure(1)
plot(vpiseq,cumseq,'-r',vpisl,cumsl,'-b');%,vpibr,cumbr,'-.k');
legend('SEQ','Streamlines');%,'Referência');
xlabel('pvi');
ylabel('cumulative oil production');
grid

saveas(gcf,sprintf('%s\\%s\\%s\\Cumulateoil.png',caso,malha,formethod));

figure(2)
plot(vpiseq,wcutseq,'-r',vpisl,wcutsl,'-b');%,vpibr,wcutbr,'-.k');
legend('SEQ','Streamlines');%,'Referência');
xlabel('pvi');
ylabel('watercut');
grid

saveas(gcf,sprintf('%s\\%s\\%s\\Watercut.png',caso,malha,formethod));

end

