function metadata = extract_meta(headers, regex, tok_no)
    % Extract a metadata header from the original metadata given a regex and
    % token number (use a regex with a token).
    if ~iscellstr(headers)
        disp('You have not provided a cell array of strings')
        headers_string = cell(size(headers));
        for i = 1:numel(headers)
            headers_string{i} = char(headers{i});
        end
        headers = headers_string;
    end
    [~, tok] = regexp(headers, regex, 'match', 'tokens');
    metadata = cell(size(tok));
    for i = 1:numel(tok)
        if numel(tok{i}) < 1
            metadata(i) = {'Missing data'};
        else
            metadata(i) = tok{i}{1}(tok_no);
        end
    end
end
