function Quu = g_Quu(in1,in2,in3,in4,in5,in6)
%G_QUU
%    QUU = G_QUU(IN1,IN2,IN3,IN4,IN5,IN6)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    25-Jun-2021 23:24:01

P3_3 = in1(11);
P3_4 = in1(15);
P4_3 = in1(12);
P4_4 = in1(16);
R1_1 = in4(1);
R1_2 = in4(3);
R2_1 = in4(2);
R2_2 = in4(4);
t2 = P3_3.*2.064394473220417e-4;
t4 = P3_4.*4.432322038968113e-3;
t5 = P4_3.*4.432322038968113e-3;
t6 = P4_4.*8.658204630614184e-3;
t3 = -t2;
t7 = -t4;
t8 = -t5;
t9 = R1_2+R2_1+t2+t6+t7+t8;
Quu = reshape([P3_4.*2.064394473220417e-4+P4_3.*2.064394473220417e-4-P4_4.*2.064394473220417e-4+R1_1.*2.0+t3,t9,t9,P3_4.*8.658204630614184e-3+P4_3.*8.658204630614184e-3-P4_4.*3.631307310595811e-1+R2_2.*2.0+t3],[2,2]);
