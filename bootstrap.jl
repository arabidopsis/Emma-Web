(pwd() != @__DIR__) && cd(@__DIR__) # allow starting app from bin/ dir

using EmmaWeb
const UserApp = EmmaWeb
EmmaWeb.main()
