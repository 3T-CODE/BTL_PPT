%3 Phuong phap tuong tu Newton
%Tac gia : Nhom 22 , Lop PPT L10

%Nhap cac ham , thu vien can thiet
pkg load symbolic;
%addpath(/home/crypt/Documents/HCMUT/HK232/PPT/Group_Exercise/Code/Code/BFGS.m);
%savepath;

%Thiet lap lai gia tri
clear all;
close all;

%Nhap gia tri ban dau cua x , ham so f(x) , va gia tri dung sai hoi tu e
x_giaTri = input("Nhap gia tri x ban dau : x = ");
global n;
n = int16(length(x_giaTri));

sym f;
x = sym('x', [1 n]);
global I ;
I = eye(n);

f = input("Nhap ham so : f(x) = ");
e = double(input("Nhap gia tri dung sai hoi tu : e = "));

%Kiem tra gia tri cua e
while e <= 0
  e = double(input("Nhap gia tri dung sai hoi tu : e = "));
end

%Chon cach tinh
fprintf("\n Ban muon tim gi ? \n");
fprintf(" \n 1/Tim cuc tieu gan dung - (Nhap 1)");
fprintf(" \n 2/Tim cuc dai gan dung - (Nhap 2)");
fprintf(" \n 3/Tim nghiem gan dung - (Nhap 3) ");
fprintf(" \n 4/Huy tinh toan - (Nhap so khac) \n");
Lua_chon_cach_tinh = input("Nhap : ") ;

%Xet lua chon da nhap
if Lua_chon_cach_tinh == 1
  %Chuyen doi ham so phu hop
  f = f ;

  %Thong bao lua chon da nhap
  fprintf("\nBan da chon tim cuc tieu \n");

elseif Lua_chon_cach_tinh == 2
  %Chuyen doi ham so phu hop
  f = -f ;

  %Thong bao lua chon da nhap
  fprintf("\nBan da chon tim cuc dai \n");

elseif Lua_chon_cach_tinh == 3
  %Chuyen doi ham so phu hop
  f = f.^2;

  %Thong bao lua chon da nhap
  fprintf("\nBan da chon tim nghiem \n");

else

end

%Lua chon phuong phap tinh
fprintf("\n Hay lua chon phuong phap tinh ");
fprintf("\n 1/ BFGS - (Nhap 1)");
fprintf("\n 2/ DFP  - (Nhap 2)");
fprintf("\n 3/ SR1 - (Nhap 3) ");
fprintf("\n 4/ Huy tinh toan  - (Nhap so khac) \n");
Lua_chon_PP_tinh = input("Nhap : ");

%Xet lua chon da nhap va khoi dong phuong phap tinh da chon
if Lua_chon_PP_tinh == 1
  %Thong bao PP da lua chon
  fprintf("\n Ban da chon phuong phap BFGS \n ");

  %Khoi dong tinh toan
  x_giaTri = BFGS( x_giaTri , f , e);

  %In gia tri ra man hinh
  fprintf("\n Dap an : \n");
  fprintf("\n X = \n");
  disp(x_giaTri);

elseif Lua_chon_PP_tinh == 2
  %Thong bao PP da lua chon
  fprintf("\n Ban da chon phuong phap DFP \n ");

  %Khoi dong tinh toan
  x_giaTri = DFP( x_giaTri , f , e);

  %In gia tri ra man hinh
  fprintf("\n Dap an : \n");
  fprintf("\n X = \n");
  disp(x_giaTri);

elseif Lua_chon_PP_tinh == 3
   %Thong bao PP da chon
   fprintf("\n Ban da chon phuong phap SR1 \n");

   %Nhap cac du lieu con thieu
   fprintf("\n Hay nhap cac gia tri can thiet sau day : \n");
   delta = input("\n Nhap ban kinh hoi tu (Delta > 0):  Delta = ");
   eta = input("\n Nhap toan tu eta thuoc (0 ; 10^-3) : eta = ");
   r = input("\n Nhap toan tu r thuoc (0;1) : r = ");

   %Kiem tra du lieu da dat dieu kien hay chua
   while 1
        if delta > 0
          if 0 < eta && eta < 1e-3
            if 0 < r && r < 1
              break;
            end
          end
        end

       %Tien hanh nhap lai du lieu
    fprintf("\n Hay nhap cac gia tri can thiet sau day : \n");
   delta = input("\n Nhap ban kinh hoi tu (Delta > 0):  Delta = ");
   eta = input("\n Nhap toan tu eta thuoc (0 ; 10^-3) : eta = ");
   r = input("\n Nhap toan tu r thuoc (0;1) : r = ");

   end

   %Tim gia tri theo SR1
   x_giaTri = SR1(x_giaTri , f , e  , delta , eta , r);

   %In gia tri da tim duoc ra man hinh
   fprintf("\n Dap an : \n");
   fprintf("\nX =  \n");
   disp(x_giaTri);

else
   %Thong bao PP da nhap
   fprintf("\n Ban da huy gia tri tinh toan \n");


end
%%%%Ket thuc %%%%










