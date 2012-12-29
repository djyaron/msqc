function updateContext(obj,ienv)

rho = obj.mcharge(ienv);
obj.charges(:,ienv+1) = rho(:);
obj.bondOrders(:,:,ienv+1) = obj.calcBO(ienv);

end

