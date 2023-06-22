## GTOC12 Verification Program
This is an automatic verification program provided for GTOC12 participants, enabling them to verify their design results more flexibly offline. We offer executable programs for Windows, Linux, and macOS operating systems. Teams with different operating systems can select the appropriate program for their system for validation.

While our testing may not cover all situations, we welcome any feedback you may have. Please feel free to contact us either by posting in the `Discussion` or by sending an email to `zhong-zh19@mails.tsinghua.edu.cn`. 

If the results of this verification program do not align with those on the website, please refer to the website version as the standard.

## Getting Started


### Structure

  `GTOC12_Asteroids_Data`: This file contains the orbit elements of all asteroids in GTOC12. (The file name cannot be changed.)

  `Result`: This is the solution file. (The file name cannot be changed.)

  `Verify_GTOC12`: This is the executable verification program. It will return result information upon successful verification, and provide the reason for failure otherwise.

  `ScoreData`: The detailed result information, which includes the number of mined asteroids, as well as the ID and the collected mass of each, will only appear after successful verification.

### Verification

* Linux & macOS

Run `./Verify_GTOC12` in the terminal, from the corresponding directory.

* Windows

To start the program, double-click on the `Verify_GTOC12.exe` file.
