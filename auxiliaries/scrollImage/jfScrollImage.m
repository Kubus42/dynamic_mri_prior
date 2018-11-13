function fid = jfScrollImage(vol, defaultSlice,fidIn)
    
    slice = max(1,floor(size(vol,3)/2));
    if nargin>1
        slice = defaultSlice;
    end
    if nargin>2
        fid = fidIn;
    else       
        fid = figure;
    end
    fid = figure(fid);
    fid.WindowScrollWheelFcn = @scrollWheelFcn;
    cRange = [min(vol(:)), max(vol(:))];
    drawSlice;

    function drawSlice
        slice = max(slice, 1);
        slice = min(slice, size(vol,3));
        subaxis(1,1,1,'s',0,'p',0,'m',0,'mt',0.075)
        imagesc(vol(:,:,slice),cRange);
        %colorbar;
        title(['Slice ', num2str(slice)]);
        colormap gray;
%         caxis([0,1]);
    end

    function scrollWheelFcn(~,evnt)
        slice = slice + evnt.VerticalScrollCount;
        %disp(num2str(slice));
        drawSlice;
    end

end