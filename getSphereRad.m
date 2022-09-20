function R = getSphereRad(Points)
% getSphereRad - get Sphere's Radius given 4 points in 3D space
% Points should be in 3 * 4 array format
    
    P1 = Points(:,1); P2 = Points(:,2); P3 = Points(:,3); P4 = Points(:,4);
    A = [(P2 - P1)';
         (P3 - P2)';
         (P4 - P3)'];
    b = 1/2 * [P2'*P2 - P1'*P1;
               P3'*P3 - P2'*P2;
               P4'*P4 - P3'*P3]; 
    Xr = A \ b;
%     disp(['Calculated Circle Center: ',num2str(Xr')])
%     disp(['Radius by P1: ',num2str(norm(Xr - P1))])
%     disp(['Radius by P2: ',num2str(norm(Xr - P2))])
%     disp(['Radius by P3: ',num2str(norm(Xr - P3))])
%     disp(['Radius by P4: ',num2str(norm(Xr - P4))])
    R = norm(Xr - P1);

%     figure(1);
%     [Xs,Ys,Zs] = sphere;
%     Xs = Xs * R + Xr(1); Ys = Ys * R + Xr(2); Zs = Zs * R + Xr(3);
%     surf(Xs,Ys,Zs); hold on; grid on;
%     plot3(Points(1,:),Points(2,:),Points(3,:),'b*'); 
%     plot3(Xr(1),Xr(2),Xr(3),'rx');
end