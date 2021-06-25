function Quu = g2_Quu(in1,in2,in3,in4,in5,in6)
%G2_QUU
%    QUU = G2_QUU(IN1,IN2,IN3,IN4,IN5,IN6)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    26-Jun-2021 00:44:54

P3_3 = in1(11);
P3_4 = in1(15);
P4_3 = in1(12);
P4_4 = in1(16);
R1_1 = in4(1);
R1_2 = in4(3);
R2_1 = in4(2);
R2_2 = in4(4);
t2 = P3_3.*2.064394473220417e-4;
t3 = P3_4.*2.064394473220417e-4;
t4 = P4_3.*2.064394473220417e-4;
t6 = P3_4.*8.658204630614184e-3;
t7 = P4_3.*8.658204630614184e-3;
t8 = P4_4.*8.658204630614184e-3;
t5 = -t2;
t9 = -t8;
Quu = reshape([P4_4.*2.064394473220417e-4+R1_1+t2-t3-t4,R2_1+t3+t5+t7+t9,R1_2+t4+t5+t6+t9,P4_4.*3.631307310595811e-1+R2_2+t2-t6-t7],[2,2]);