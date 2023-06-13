void geomet(double hx, double hy, double& rr, double& theta);

void geomet(double hx, double hy, double& rr, double& theta)
{
    rr = sqrt(hx*hx + hy*hy);
    theta = atan2(hx,hy);
    
}


