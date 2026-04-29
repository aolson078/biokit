function compileReport() {
    window.location.href = '/compile_report';
}

function displayReport(reportId) {
    window.location.href = '/display_report/' + reportId;
}

async function postJson(url, payload) {
    const response = await fetch(url, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload),
    });

    if (!response.ok) {
        throw new Error(`Request failed with status: ${response.status}`);
    }

    return response.json();
}

function displayResults(results) {
    const searchResultsContainer = document.getElementById('searchResults');
    if (!searchResultsContainer) return;

    searchResultsContainer.innerHTML = '';

    results.forEach(function(result) {
        const card = document.createElement('div');
        card.className = 'card';

        const title = document.createElement('h3');
        title.textContent = result.title;
        card.appendChild(title);

        let sequence = result.description;
        const characterLimit = 500;
        if (sequence.length > characterLimit) {
            sequence = sequence.substring(0, characterLimit) + '...';
        }

        const sequenceElement = document.createElement('p');
        sequenceElement.textContent = sequence;
        card.appendChild(sequenceElement);

        const selectButton = document.createElement('button');
        selectButton.className = 'select-button';
        selectButton.textContent = 'Select';
        selectButton.addEventListener('click', function() {
            const selectedResult = result.title + ' ' + result.description;
            sendSelectedResult(selectedResult);
        });
        card.appendChild(selectButton);

        searchResultsContainer.appendChild(card);
    });
}

function mountReactSearchResults(initialResults) {
    const appRoot = document.getElementById('searchResultsApp');
    if (!appRoot || !window.React || !window.ReactDOM) {
        displayResults(initialResults);
        return;
    }

    const { createElement: h, useState } = window.React;

    function SearchResultsApp() {
        const [results, setResults] = useState(initialResults || []);

        window.__setSearchResults = setResults;

        return h(
            'div',
            { className: 'card-container' },
            results.map((result, index) => {
                const description = result.description.length > 500
                    ? `${result.description.substring(0, 500)}...`
                    : result.description;

                return h('div', { className: 'card', key: `${result.title}-${index}` }, [
                    h('h3', { key: 'title' }, result.title),
                    h('p', { key: 'description' }, description),
                    h(
                        'button',
                        {
                            key: 'select',
                            className: 'select-button',
                            onClick: function() {
                                sendSelectedResult(`${result.title} ${result.description}`);
                            }
                        },
                        'Select'
                    )
                ]);
            })
        );
    }

    window.ReactDOM.createRoot(appRoot).render(h(SearchResultsApp));
}

async function searchReports() {
    const searchInput = document.getElementsByName('search')[0];
    const statusEl = document.getElementById('search-status');
    const searchQuery = searchInput ? searchInput.value.trim() : '';

    if (searchQuery.length < 3) {
        if (statusEl) statusEl.textContent = 'Enter at least 3 characters to search.';
        return;
    }

    if (statusEl) statusEl.textContent = 'Searching NCBI records...';

    try {
        const results = await postJson('/search', { search: searchQuery });
        if (statusEl) statusEl.textContent = `Found ${results.length} result(s).`;
        if (window.__setSearchResults) {
            window.__setSearchResults(results);
        } else {
            displayResults(results);
        }
    } catch (error) {
        console.error('Search request failed:', error);
        if (statusEl) statusEl.textContent = 'Search failed. Please try again.';
    }
}

async function sendSelectedResult(selectedResult) {
    try {
        await postJson('/employee.html/', { selected_result: selectedResult });
        window.location.reload();
    } catch (error) {
        console.error('Failed to send selected result:', error);
    }
}

document.addEventListener('DOMContentLoaded', function() {
    mountReactSearchResults([]);

    const userForm = document.getElementById('select-user-form');
    if (userForm) {
        userForm.addEventListener('submit', function(event) {
            event.preventDefault();
            const dropdown = document.getElementById('user-dropdown');
            const userId = dropdown.value;
            const username = dropdown.options[dropdown.selectedIndex].text;
            sendSelectedUser(userId, username);
        });
    }


    const clearSearchBtn = document.getElementById('clearSearchBtn');
    if (clearSearchBtn) {
        clearSearchBtn.addEventListener('click', function() {
            const searchInput = document.getElementsByName('search')[0];
            const statusEl = document.getElementById('search-status');
            if (searchInput) searchInput.value = '';
            if (window.__setSearchResults) {
                window.__setSearchResults([]);
            } else {
                displayResults([]);
            }
            if (statusEl) statusEl.textContent = 'Search cleared.';
        });
    }

    const deleteButtons = document.querySelectorAll('.delete-btn');
    deleteButtons.forEach(function(button) {
        button.addEventListener('click', async function() {
            const userId = button.dataset.userid;
            try {
                const response = await fetch(`/delete_user/${userId}`, { method: 'DELETE' });
                if (response.ok) {
                    button.closest('tr').remove();
                    alert('User deleted successfully');
                } else {
                    const data = await response.json();
                    alert(`Error: ${data.error}`);
                }
            } catch (error) {
                console.error('An error occurred while deleting user:', error);
                alert('An error occurred while deleting user');
            }
        });
    });
});

function sendSelectedUser(userId, username) {
    const xhr = new XMLHttpRequest();
    xhr.onreadystatechange = function() {
        if (xhr.readyState === XMLHttpRequest.DONE) {
            if (xhr.status === 200) {
                document.querySelector('.user-id').textContent = userId;
            } else {
                console.error('Failed to send selected user with status:', xhr.status);
            }
        }
    };
    xhr.open('POST', '/process_selected_user', true);
    xhr.setRequestHeader('Content-Type', 'application/json');
    xhr.send(JSON.stringify({ userId: userId, username: username }));
}

function logout() {
    var xhr = new XMLHttpRequest();
    xhr.onreadystatechange = function() {
        if (xhr.readyState === XMLHttpRequest.DONE) {
            if (xhr.status === 200) {
                window.location.href = '/index';
            } else {
                console.error('Failed to logout with status:', xhr.status);
            }
        }
    };
    xhr.open('GET', '/logout', true);
    xhr.send();
}

function showUsernameAndPasswordFields() {
    const usernameInput = document.getElementById('usernameInput');
    const passwordInput = document.getElementById('passwordInput');
    const usernameLabel = document.querySelector("label[for='usernameInput']");
    const passwordLabel = document.querySelector("label[for='passwordInput']");
    const selectedUser = document.querySelector('.user-id').textContent;

    if (selectedUser !== 'None selected') {
        usernameInput.style.display = 'block';
        passwordInput.style.display = 'block';
    } else {
        usernameInput.style.display = 'none';
        passwordInput.style.display = 'none';
        usernameLabel.style.display = 'none';
        passwordLabel.style.display = 'none';
    }
}
