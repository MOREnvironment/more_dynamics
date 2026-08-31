"""Keep pytest and ROS launch-test discovery out of runtime plugin sources."""

collect_ignore = [
    "more_dynamics/plugins",
    "more_dynamics/plugin_types",
]

collect_ignore_glob = [
    "more_dynamics/plugins/**",
    "more_dynamics/plugin_types/**",
]
