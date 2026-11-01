function deleteWaitbar(wbar)

    try

        if ishghandle(wbar)
            delete(wbar)
        end

    catch
    end

end
