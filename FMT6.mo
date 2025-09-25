

model FullMechanismT6
  parameter Real l2 = 45.0 "mm";
  parameter Real l3 = 110.0 "mm";
  parameter Real l4 = 59.0 "mm";
  parameter Real l5 = 58.0 "mm";
  parameter Real l6 = 46.0 "mm";

  parameter Real omega = -0.3 * 3.141592653589793 "rad/s";

  Real theta1(start=0);
  Real d1(start=103.2389, fixed=false);
  Real theta2(start=1.4652, fixed=false);
  Real theta3(start=0.1663, fixed=false);
  Real d7(start=-33.7662, fixed=false);
  

  Real d1_dot = der(d1);
  Real theta2_dot = der(theta2);
  Real theta3_dot = der(theta3);
  Real d7_dot = der(d7);

  Real d1_ddot = der(d1_dot);
  Real theta2_ddot = der(theta2_dot);
  Real theta3_ddot = der(theta3_dot);
  Real d7_ddot = der(d7_dot);

  Real xP, yP, vxP, vyP, axP, ayP, aP;

equation
  theta1 = omega*time;

  l2*cos(theta1) - d1*cos(theta2) = 0;
  l4 + l2*sin(theta1) - d1*sin(theta2) = 0;
  l2*cos(theta1) + (l3 - d1)*cos(theta2) - l6*cos(theta3) - d7 = 0;
  l2*sin(theta1) + (l3 - d1)*sin(theta2) + l6*sin(theta3) - l5 = 0;

  xP = l3*cos(theta2);
  yP = l3*sin(theta2);

  vxP = der(xP);
  vyP = der(yP);
  axP = der(vxP);
  ayP = der(vyP);
  
  aP = sqrt(axP^2 + ayP^2);
  
end FullMechanismT6;
