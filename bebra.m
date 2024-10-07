for n = 1:16
    subplot(4, 4, n)
    ord = n
    m = magic(ord)
    imagesc(m)
    title(num2str(ord))
    axis equal
    axis off
end