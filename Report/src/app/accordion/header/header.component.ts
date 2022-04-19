import {Component, EventEmitter, Input, OnInit, Output} from '@angular/core';

@Component({
  selector: 'accordion-header',
  templateUrl: './header.component.html',
  styleUrls: ['./header.component.css']
})
export class HeaderComponent implements OnInit {

  @Input()
  title: string = "undefined";

  @Input()
  disabled: boolean = false;

  @Input()
  sub: boolean = false;

  @Input()
  information: string | undefined;

  @Output()
  togglePanel = new EventEmitter<boolean>();

  informationVisible: boolean = false;
  panelVisible: boolean = false;

  constructor() {
  }

  ngOnInit(): void {
  }

  onPanelClick() {
    this.panelVisible = !this.panelVisible;

    this.togglePanel.emit(this.panelVisible);
  }
}
