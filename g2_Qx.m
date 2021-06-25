function Qx = g2_Qx(in1,in2,in3,in4,in5,in6)
%G2_QX
%    QX = G2_QX(IN1,IN2,IN3,IN4,IN5,IN6)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    26-Jun-2021 00:44:54

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
p1 = in2(:,1);
p2 = in2(:,2);
p3 = in2(:,3);
p4 = in2(:,4);
x1 = in5(1,:);
x2 = in5(2,:);
x3 = in5(3,:);
x4 = in5(4,:);
t2 = cos(x1);
Qx = [p1+Q1_1.*x1+Q2_1.*x2+Q3_1.*x3+Q4_1.*x4-p3.*t2.*2.956624792929091e-2+p4.*t2.*2.956624792929091e-2,p2+Q1_2.*x1+Q2_2.*x2+Q3_2.*x3+Q4_2.*x4,p1./1.0e+3+p3+Q1_3.*x1+Q2_3.*x2+Q3_3.*x3+Q4_3.*x4,p2./1.0e+3+p3.*2.155200121739496e-6+p4.*9.999096095057606e-1+Q1_4.*x1+Q2_4.*x2+Q3_4.*x3+Q4_4.*x4];