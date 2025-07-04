    function compileReport() {
        window.location.href = '/compile_report';
    }

    function displayReport(reportId) {
        window.location.href = '/display_report/' + reportId;
    }


    function searchReports() {
        var searchQuery = document.getElementsByName("search")[0].value;
        var xhr = new XMLHttpRequest();
        xhr.onreadystatechange = function() {
            if (xhr.readyState === XMLHttpRequest.DONE) {
                if (xhr.status === 200) {
                    var searchResultsContainer = document.getElementById("searchResults");
                    searchResultsContainer.innerHTML = "";

                    var results = JSON.parse(xhr.responseText);
                    displayResults(results);
                } else {
                    console.error("Search request failed with status:", xhr.status);
                }
            }
        };
        xhr.open("POST", "/search", true);
        xhr.setRequestHeader("Content-Type", "application/json");
        xhr.send(JSON.stringify({ search: searchQuery }));
    }

    function displayResults(results) {
        var searchResultsContainer = document.getElementById("searchResults");
        searchResultsContainer.innerHTML = "";

        console.log(results)

        results.forEach(function(result) {
            var card = document.createElement("div");
            card.className = "card";

            var title = document.createElement("h3");
            title.textContent = result.title;
            card.appendChild(title);

            var sequence = result.description;

            var characterLimit = 500;
            if (sequence.length > characterLimit) {
                sequence = sequence.substring(0, characterLimit) + "...";
            }

            var sequenceElement = document.createElement("p");
            sequenceElement.textContent = sequence;
            card.appendChild(sequenceElement);

            var selectButton = document.createElement("button");
            selectButton.className = "select-button";
            selectButton.textContent = "Select";
            selectButton.addEventListener("click", function() {
                var selectedResult = result.title + " " + result.description;
                sendSelectedResult(selectedResult);
            });
            card.appendChild(selectButton);

            searchResultsContainer.appendChild(card);
        });
    }

function sendSelectedResult(selectedResult) {
    var xhr = new XMLHttpRequest();
    xhr.onreadystatechange = function() {
        if (xhr.readyState === XMLHttpRequest.DONE) {
            if (xhr.status === 200) {
                console.log("Selected result sent successfully");
                // Reload the page to reflect the updated data
                window.location.reload();
            } else {
                console.error("Failed to send selected result with status:", xhr.status);
            }
        }
    };
    xhr.open("POST", "/employee.html/", true);
    xhr.setRequestHeader("Content-Type", "application/json");
    xhr.send(JSON.stringify({ selected_result: selectedResult }));
}

document.addEventListener("DOMContentLoaded", function() {
    document.getElementById("select-user-form").addEventListener("submit", function(event) {
        event.preventDefault(); // Prevent the default form submission behavior
        // Get the selected user from the dropdown menu
        var userId = document.getElementById("user-dropdown").value;
        var username = document.getElementById("user-dropdown").options[document.getElementById("user-dropdown").selectedIndex].text;
        // Send the selected user information to the server
        sendSelectedUser(userId, username);
    });
});


function sendSelectedUser(userId, username) {
    var xhr = new XMLHttpRequest();
    xhr.onreadystatechange = function() {
        if (xhr.readyState === XMLHttpRequest.DONE) {
            if (xhr.status === 200) {
                console.log("Selected user sent successfully");
                document.querySelector(".user-id").textContent = userId;
                //window.location.reload();
            } else {
                console.error("Failed to send selected user with status:", xhr.status);
            }
        }
    };
    xhr.open("POST", "/process_selected_user", true);
    xhr.setRequestHeader("Content-Type", "application/json");
    xhr.send(JSON.stringify({ userId: userId, username: username }));
}

// added
document.addEventListener('DOMContentLoaded', function() {
    const deleteButtons = document.querySelectorAll('.delete-btn');
    deleteButtons.forEach(button => {
        button.addEventListener('click', function() {
            const userId = button.dataset.userid;
            fetch(`/delete_user/${userId}`, {
                method: 'DELETE',
            })
            .then(response => {
                if (response.ok) {
                    // Optionally, remove the deleted user from the UI
                    button.closest('tr').remove();
                    alert('User deleted successfully');
                } else {
                    response.json().then(data => {
                        alert(`Error: ${data.error}`);
                    }).catch(error => {
                        console.error('Error parsing JSON:', error);
                        alert('An error occurred while deleting user');
                    });
                }
            })
            .catch(error => {
                console.error('Error:', error);
                alert('An error occurred while deleting user');
            });
        });
    });
});

function logout() {
    var xhr = new XMLHttpRequest();
    xhr.onreadystatechange = function() {
        if (xhr.readyState === XMLHttpRequest.DONE) {
            if (xhr.status === 200) {
                console.log("Logged out successfully");
                window.location.href = '/index'; // Redirect to login page after logout
            } else {
                console.error("Failed to logout with status:", xhr.status);
            }
        }
    };
    xhr.open("GET", "/logout", true);
    xhr.send();
}

function showUsernameAndPasswordFields() {
            const usernameInput = document.getElementById("usernameInput");
            const passwordInput = document.getElementById("passwordInput");
            var usernameLabel = document.querySelector("label[for='usernameInput']");
            var passwordLabel = document.querySelector("label[for='passwordInput']");
            const selectedUser = document.querySelector(".user-id").textContent;

            if (selectedUser !== "None selected") {
                usernameInput.style.display = "block";
                passwordInput.style.display = "block";

            } else {
                usernameInput.style.display = "none";
                passwordInput.style.display = "none";
                usernameLabel.style.display = "none";
                passwordLabel.style.display = "none";
            }
        }