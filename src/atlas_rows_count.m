function n = atlas_rows_count(rows)
%RS3_ROWS_COUNT  Number of rows (handles packed struct or double matrix).

if isstruct(rows)
    n = double(rows.n);
elseif isnumeric(rows)
    n = size(rows, 1);
else
    n = 0;
end
end
