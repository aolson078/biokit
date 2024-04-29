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
                window.location.reload();
            } else {
                console.error("Failed to send selected user with status:", xhr.status);
            }
        }
    };
    xhr.open("POST", "/process_selected_user", true);
    xhr.setRequestHeader("Content-Type", "application/json");
    xhr.send(JSON.stringify({ userId: userId, username: username }));
}