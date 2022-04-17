import { ComponentFixture, TestBed } from '@angular/core/testing';

import { TfGroupsContainerComponent } from './tf-groups-container.component';

describe('TfGroupsContainerComponent', () => {
  let component: TfGroupsContainerComponent;
  let fixture: ComponentFixture<TfGroupsContainerComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ TfGroupsContainerComponent ]
    })
    .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(TfGroupsContainerComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
