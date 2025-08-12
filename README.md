# DNAAT: A Lightweight DNA Alignment Tool

DNAAT is a lightweight and simple tool for performing local and global DNA sequence alignments. It is written in Python and uses C-implemented functions for improved performance. This tool is designed for quick and easy use without the need for heavy dependencies or complex installations.

-----

## Features

  * **Global DNA sequence alignments**: Compare sequences of similar length.
  * **Local DNA sequence alignments**: Find the best region of similarity between two sequences.
  * **Lightweight and easy to use**: Designed for simplicity and minimal setup.
  * **Performance Testing**: Includes a script to compare the performance of the Python and C implementations.

-----

## Getting Started

To get started with DNAAT, follow these steps:

1.  **Clone the repository:**

    ```bash
    git clone https://github.com/no1berg/DNAAT.git
    cd DNAAT
    ```

2.  **Run an alignment:**

    Use the `dnaat.py` script with either the `global` or `local` command.

      * For **global alignment**:
        ```bash
        python dnaat.py global
        ```
      * For **local alignment**:
        ```bash
        python dnaat.py local
        ```
    *You can also provide a path to a FASTA file as an argument.*

4.  **Test Performance:**
    To see the performance boost from the C implementation, run:

    ```bash
    python test_performance.py
    ```

-----

## Learn More

For more details, check out the project's GitHub Pages site:
[https://no1berg.github.io/dnaat/](https://no1berg.github.io/dnaat/)

-----

## Contact

If you have any questions or suggestions, feel free to reach out:

  * **Email**: [no1sebastian@gmail.com](mailto:no1sebastian@gmail.com)
  * **GitHub**: [https://github.com/no1berg/DNAAT](https://github.com/no1berg/DNAAT)
