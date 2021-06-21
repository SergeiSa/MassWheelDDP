function Qxx = g_Qxx(in1,in2,in3,R1,in5,u1)
%G_QXX
%    QXX = G_QXX(IN1,IN2,IN3,R1,IN5,U1)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    21-Jun-2021 13:27:03

P1_1 = in1(1);
P1_2 = in1(5);
P1_3 = in1(9);
P1_4 = in1(13);
P2_1 = in1(2);
P2_2 = in1(6);
P2_3 = in1(10);
P2_4 = in1(14);
P3_1 = in1(3);
P3_2 = in1(7);
P3_3 = in1(11);
P3_4 = in1(15);
P4_1 = in1(4);
P4_2 = in1(8);
P4_3 = in1(12);
P4_4 = in1(16);
Q1_1 = in3(1);
Q1_2 = in3(5);
Q1_3 = in3(9);
Q1_4 = in3(13);
Q2_1 = in3(2);
Q2_2 = in3(6);
Q2_3 = in3(10);
Q2_4 = in3(14);
Q3_1 = in3(3);
Q3_2 = in3(7);
Q3_3 = in3(11);
Q3_4 = in3(15);
Q4_1 = in3(4);
Q4_2 = in3(8);
Q4_3 = in3(12);
Q4_4 = in3(16);
p3 = in2(:,3);
p4 = in2(:,4);
x1 = in5(1,:);
x2 = in5(2,:);
x3 = in5(3,:);
x4 = in5(4,:);
t2 = cos(x1);
t3 = sin(x1);
t4 = -x3;
t5 = P1_1./5.0e+2;
t6 = P2_2./5.0e+2;
t7 = P1_2./1.0e+3;
t8 = P1_3./1.0e+3;
t9 = P1_4./1.0e+3;
t10 = P2_1./1.0e+3;
t11 = P2_3./1.0e+3;
t12 = P2_4./1.0e+3;
t13 = P3_2./1.0e+3;
t14 = P4_1./1.0e+3;
t15 = x3./1.0e+3;
t16 = x4./1.0e+3;
t23 = P1_2./1.0e+6;
t24 = P2_1./1.0e+6;
t26 = u1.*8.764480495073951e-4;
t47 = u1.*3.675880099068387e-2;
t17 = P3_3+t8;
t18 = P3_4+t9;
t19 = P4_3+t11;
t20 = P4_4+t12;
t21 = t15+x1;
t22 = t16+x2;
t25 = P2_4+P4_2+Q2_4+Q4_2+t6;
t27 = t3.*2.956624792929091e-2;
t28 = P2_3+P3_2+Q2_3+Q3_2+t7+t10;
t29 = P2_3.*t2.*2.956624792929091e-2;
t30 = P2_4.*t2.*2.956624792929091e-2;
t31 = P3_2.*t2.*2.956624792929091e-2;
t32 = P3_3.*t2.*2.956624792929091e-2;
t33 = P3_4.*t2.*2.956624792929091e-2;
t34 = P4_2.*t2.*2.956624792929091e-2;
t35 = P4_3.*t2.*2.956624792929091e-2;
t36 = P4_4.*t2.*2.956624792929091e-2;
t41 = P3_1.*t2.*2.956624792929091e-5;
t42 = P3_2.*t2.*2.956624792929091e-5;
t43 = P4_1.*t2.*2.956624792929091e-5;
t44 = P4_2.*t2.*2.956624792929091e-5;
t37 = -t29;
t38 = -t31;
t39 = -t32;
t40 = -t33;
t45 = -t41;
t46 = -t42;
t48 = t2.*t17.*2.956624792929091e-2;
t49 = t2.*t18.*2.956624792929091e-2;
t50 = t2.*t19.*2.956624792929091e-2;
t51 = t2.*t20.*2.956624792929091e-2;
t54 = t4+t26+t27-2.155200121739496e-6;
t55 = Q3_4+Q4_3+t13+t14+t18+t19+t23+t24;
t56 = t27+t47+x4-9.039049423938655e-5;
t52 = -t48;
t53 = -t50;
t57 = P1_2+P2_1+Q1_2+Q2_1+t30+t34+t37+t38;
t58 = P1_3+P3_1+Q1_3+Q3_1+t5+t35+t39+t43+t45+t49+t52;
t59 = P1_4+P4_1+Q1_4+Q4_1+t7+t10+t36+t40+t44+t46+t51+t53;
Qxx = reshape([P1_1.*2.0+Q1_1.*2.0-P3_1.*t2.*5.913249585858181e-2+P4_1.*t2.*5.913249585858181e-2-p4.*t3.*2.956624792929091e-2+p3.*t27-t2.*(P1_3+t35+t39).*5.913249585858181e-2+t2.*(P1_4+t36+t40).*5.913249585858181e-2-t3.*(P1_4.*t21+P2_4.*t22-P3_4.*t54+P4_4.*t56).*2.956624792929091e-2+t27.*(P1_3.*t21+P2_3.*t22-P3_3.*t54+P4_3.*t56)+t3.*t21.*(P3_1.*2.956624792929091e-2-P4_1.*2.956624792929091e-2)+t3.*t22.*(P3_2.*2.956624792929091e-2-P4_2.*2.956624792929091e-2)-t3.*t54.*(P3_3.*2.956624792929091e-2-P4_3.*2.956624792929091e-2)+t3.*t56.*(P3_4.*2.956624792929091e-2-P4_4.*2.956624792929091e-2),t57,t58,t59,t57,P2_2.*2.0+Q2_2.*2.0,t28,t25,t58,t28,P1_1./5.0e+5+P1_3./5.0e+2+P3_1./5.0e+2+P3_3.*2.0+Q3_3.*2.0,t55,t59,t25,t55,P2_2./5.0e+5+P2_4./5.0e+2+P4_2./5.0e+2+P4_4.*2.0+Q4_4.*2.0],[4,4]);
