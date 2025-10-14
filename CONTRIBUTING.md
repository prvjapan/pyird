# How to Contribute to PyIRD

Thank you for considering contributing to `PyIRD`! To ensure smooth and efficient development, please follow these guidelines.

---

## 1. Branching Strategy

We follow the branch flow outlined below:

- **`master` branch**: The latest stable version. Do not commit directly to this branch.  
- **`release` branch**: Code that is ready for release.  
- **`develop` branch**: The main development branch where all new changes are merged.  
- **`contributor-defined-branch`**: A feature or bugfix branch created by each contributor.

### Workflow:
1. **Fork** the repository and clone it locally.  
2. Start a new branch from `develop`. Use descriptive branch names like `feature/<name>` or `bugfix/<name>`.  

   ```bash
   git checkout develop
   git checkout -b feature/<name>
   ```

3. Add your code, **write unit tests**, and commit your changes.  
4. Push your branch to the remote repository and open a **Pull Request (PR)**.  

   **Base branch**: `develop`.

   ```bash
   git push origin feature/<name>
   ```

5. Your PR will be reviewed. After approval, it will be merged into `develop`.  
6. Once stable, the code will be merged into `release`, and finally into `master` for deployment.

---

## 2. Creating Issues

If you encounter a bug or have suggestions, feel free to open an **Issue**.  

- Clearly describe the problem (include reproduction steps, expected vs. actual behavior).  
- Use a concise and descriptive title.  
- Add appropriate labels (e.g., `bug`, `enhancement`, `question`, etc.).

---

## 3. Testing Guidelines

All contributions must include tests to ensure stability and functionality:

- Write **unit tests** for any new code or functionality.  
- Verify that all tests pass locally before submitting your Pull Request.  
- Use a testing framework of `pytest`.

See also [docs](https://secondearths.sakura.ne.jp/pyird/developers/pytest.html).

### Example: Running tests with `pytest`:

```bash
pytest tests/unittests
pytest tests/integration
```

---

## 4. Code Style

Maintain clean and readable code. Follow these guidelines:

- Adhere to **PEP8** (while not all existing code strictly adheres to it at this time, we aim to improve consistency going forward).
- Use linters (e.g., `flake8`) and formatters (e.g., `black`) to ensure consistency.

### Example: Formatting code with `black`:

```bash
black .
flake8 .
```

---

## 5. Pull Requests (PRs)

When creating a Pull Request:

- Set the base branch to `develop`.  
- Provide a clear description of the changes.  
- Link any related Issues (e.g., `Fixes #123`).  
- Ensure that all tests pass and include new tests if applicable.

---

## Questions or Feedback

If you have any questions or need clarification, please open an Issue or leave a comment in your Pull Request.

---

Thank you for contributing! We appreciate your support.  
Happy Coding 🎉  

---