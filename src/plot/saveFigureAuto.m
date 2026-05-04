function saveFigureAuto(figHandle, saveDir, fileName)
% saveFigureAuto - Saves a figure in both PDF and PNG formats with automatic handling of directories and compatibility.
% 
% SYNTAX:
%  saveFigureAuto(figHandle, saveDir, fileName)
%
% INPUT:
%  figHandle: Handle to the figure to be saved
%  saveDir: Directory where the figure should be saved
%  fileName: Name of the file (without extension)

    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    pdfFile = fullfile(saveDir, [fileName '.pdf']);
    pngFile = fullfile(saveDir, [fileName '.png']);

    try
        % PDF vectorial
        exportgraphics(figHandle, pdfFile, 'ContentType', 'vector');
        
        % PNG opcional por si también lo quieres
        exportgraphics(figHandle, pngFile, 'Resolution', 300);

    catch
        % Fallback para MATLAB antiguos
        print(figHandle, pdfFile, '-dpdf', '-bestfit');
        print(figHandle, pngFile, '-dpng', '-r300');
    end

end