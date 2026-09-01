function [Vout,Vin,Vratio,V_SUM,f] = GetData_out_in_ratio_f_VSUM(fileName)
 fname = sprintf('%s.txt',fileName);
 file_exp = fopen(fname,'r');
 Data = fscanf(file_exp,'%f %f %f %f', [4,100]);
 Vin = Data(1,:)';
 Vout = Data(2,:)';
 f = Data(3,:)';
 V_SUM = Data(4,:)';
 Vratio= -Vin./Vout;
 fclose(file_exp);
  
